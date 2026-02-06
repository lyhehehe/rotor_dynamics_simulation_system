%% calculateResponse - Compute time-domain response of rotor-bearing systems
%
% This function calculates the dynamic response of rotor-bearing systems 
% using various numerical integration methods, supporting both custom and 
% standard rotational speed profiles.
%
%% Syntax
%  [q, dq, t] = calculateResponse(Parameter, tSpan, samplingFrequency)
%  [q, dq, t, convergenceStr] = calculateResponse(Parameter, tSpan, samplingFrequency, q0, NameValueArgs)
%
%% Description
% |calculateResponse| computes the time-domain response of rotor systems 
% under specified operating conditions. The function:
% * Supports multiple numerical integration methods
% * Handles custom and standard rotational speed profiles
% * Supports Mass Matrix and Jacobian Matrix for ODE solvers
% * Provides automatic initial condition generation
% * Includes data reduction for large simulations
% * Monitors numerical convergence during integration
%
%% Input Arguments
% * |Parameter| - System configuration structure containing:
%   * |Mesh|: Discretization data with |dofNum| (total DOF count)
%   * |Status|: Operational parameters including:
%     * |vmax|, |vmin|: Max/min rotational speeds [rad/s]
%     * |acceleration|: Rotational acceleration [rad/s²]
%     * |isUseCustomize|: Custom speed profile flag
%   * |Shaft|: Shaft configuration data
% * |tSpan| - Simulation time range [start, end] [s]
% * |samplingFrequency| - Output sampling rate [Hz]
% * |q0| - Initial displacement vector (optional)
% * |NameValueArgs| - Optional name-value pairs:
%   * |isPlotStatus|: Plot rotational status (default: true)
%   * |reduceInterval|: Data downsampling factor (default: 1)
%   * |calculateMethod|: Solver type ('RK', 'ode45', 'ode15s', 'ode23s') (default: 'RK')
%   * |options|: ODE solver options (for MATLAB ODE solvers)
%   * |isUseMassMatrix|: Enable mass matrix for ODE solvers (default: false)
%   * |isUseJacobian|: Enable Jacobian matrix for ODE solvers (default: false)
%   * |customJacobianFun|: Handle for custom force Jacobian calculation
%   * |isUseBalanceAsInitial|: Use static balance position as initial condition (default: false)
%   * |isFreshInitial|: Force recalculation of balance position (default: false)
%
%% Output Arguments
% * |q| - Displacement time history [n×m matrix, DOFs × time]
% * |dq| - Velocity time history [n×m matrix]
% * |t| - Time vector [1×m]
% * |convergenceStr| - Convergence status message (empty if converged)
%
%% Numerical Integration Methods
% 1. ​**Runge-Kutta (RK)​**:
%    * Classic 4th-order explicit method
%    * Step-by-step time integration
%    * Real-time convergence monitoring
% 2. ​**MATLAB ODE Solvers**:
%    * |ode45|: Non-stiff problems
%    * |ode15s|: Stiff problems
%    * |ode23s|: Moderately stiff problems
%    * Supports optional Mass Matrix and Jacobian for improved performance
%
%% Initial Condition Handling
% * Automatic static balance calculation when |isUseBalanceAsInitial=true|
% * Balance position caching in 'balancePosition.mat'
% * Zero initial conditions by default
% * Custom initial vectors supported
%
%% Data Management
% * Time vector generation with specified sampling rate
% * Optional downsampling via |reduceInterval|
% * Convergence warnings for non-convergent solutions
%
%% Rotational Status Processing
% When |isPlotStatus=true|:
% 1. Computes rotational parameters for each shaft:
%    * Angular position, velocity, acceleration
% 2. Supports both standard and custom speed profiles
% 3. Generates rotational status plots
%
%% Examples
% % Basic simulation with Runge-Kutta (After modeling)
% [q, dq, t] = calculateResponse(sysParams, [0 10], 1000);
%
% % ODE15s with Mass Matrix and Jacobian (After modeling)
% [q, dq, t] = calculateResponse(sysParams, [0 5], 2000, ...
%     'calculateMethod', 'ode15s', 'isUseMassMatrix', true, 'isUseJacobian', true);
%
% % Use static balance as initial condition (After modeling)
% [q, dq, t] = calculateResponse(sysParams, [0 3], 1500, ...
%     'isUseBalanceAsInitial', true);
%
%% Dependencies
% * |dynamicEquation| - System equation of motion
% * |rungeKutta| - Custom RK4 implementation
% * |calculateBalance| - Static equilibrium calculation
% * |plotRunningStatus| - Rotational parameter visualization
%
%% See Also
% dynamicEquation, rungeKutta, ode45, ode15s, calculateBalance
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function [q, dq, t, convergenceStr] = calculateResponse(Parameter, tSpan, samplingFrequency, q0, NameValueArgs)
% check input
arguments % name value pair
    Parameter
    tSpan
    samplingFrequency
    q0 = 0; % initial value of the displacement
    NameValueArgs.isPlotStatus = true; 
    NameValueArgs.reduceInterval = 1;
    NameValueArgs.calculateMethod = 'RK';
    NameValueArgs.options = odeset();
    NameValueArgs.isUseMassMatrix = false;
    NameValueArgs.isUseJacobian = false;
    NameValueArgs.customJacobianFun = @(yn, dyn, tn, omega, domega, ddomega, Parameter) customForceJacobian(yn, dyn, tn, omega, domega, ddomega, Parameter)
    NameValueArgs.isUseBalanceAsInitial = false;
    NameValueArgs.isFreshInitial = false
end


%%

% initial value
if q0 ==0
    isUseBalanceAsInitial = NameValueArgs.isUseBalanceAsInitial;
    isFreshInitial = NameValueArgs.isFreshInitial;
    if isUseBalanceAsInitial && isFreshInitial
        q0 = calculateBalance(Parameter);
    elseif isUseBalanceAsInitial && ~isFreshInitial
        % detect exist initial value
        isExistInitialValue = exist('balancePosition.mat', 'file');
        if isExistInitialValue
            % check the dimension of the initial value
            load('balancePosition.mat','qBalance')
            if length(qBalance)== Parameter.Mesh.dofNum
                q0 = qBalance;
            else
                q0 = calculateBalance(Parameter);
            end % end if length(qBalance)~=
        else
            q0 = calculateBalance(Parameter);
        end % end if isExistInitialValue
    else
        q0 = zeros(Parameter.Mesh.dofNum,1);
    end % end if isUseBalanceAsInitial && isFreshInitial
end % end if q0 ~= 0

%%

% generate time series
tStart = tSpan(1);
tEnd = tSpan(2);
tNum = floor((tEnd - tStart) * samplingFrequency);
step = 1/samplingFrequency;
t = linspace(tStart,tEnd,tNum);

%%

% initial the response
dofNum   = Parameter.Mesh.dofNum;
dyn = zeros(dofNum,1);
if q0==0
    yn = zeros(dofNum,1); % initial condition of differential euqtion
else
    if length(q0)~=dofNum
        error('Wrong dimension of the initial value. The dimension should equal to total dof number.')
    else
        yn = q0;
    end % end if 
end % end if


% update the ode options
ode_opt = NameValueArgs.options;

% enable mass matrix function in ode() solver
if NameValueArgs.isUseMassMatrix
    % Original system is second-order: M * q'' + C*q' + K*q = f
    % ODE solvers expect first-order system M_full * y' = f_full,
    % where y = [q; q'], so construct mass matrix for first-order form:
    % M_full = [I, 0; 0, M]
    dof = Parameter.Mesh.dofNum;
    M = Parameter.Matrix.mass; % assume sparse mass matrix (dof x dof)
    % Ensure M is sparse and correct size
    if ~issparse(M)
        M = sparse(M);
    end
    I = speye(dof);
    % construct block-diagonal mass matrix (sparse)
    mass_full = blkdiag(I, M);
    % set odeset Mass option
    mass_opt = odeset('Mass', mass_full);
    ode_opt = odeset(mass_opt, ode_opt);
end

% enable Jacobian matrix in ode() solver
if NameValueArgs.isUseJacobian
    % save jacobian matrix of customize force in Parameter
    if Parameter.ComponentSwitch.hasCustom
        Parameter.Custom.jacobian = NameValueArgs.customJacobianFun;
    end
    % construct full jacobian matrix
    if NameValueArgs.isUseMassMatrix
        jacobian_full_enable_mass_function = @(tn, yn) jacobian_enable_mass_function(tn, yn(1:dofNum), yn(dofNum+1:end), Parameter);
        jacobian_opt = odeset('Jacobian', jacobian_full_enable_mass_function);
    else
        jacobian_full_disable_mass_function = @(tn, yn) jacobian_disable_mass_function(tn, yn(1:dofNum), yn(dofNum+1:end), Parameter);
        jacobian_opt = odeset('Jacobian', jacobian_full_disable_mass_function);
    end
    % set odeset jacobian option
    
    ode_opt = odeset(jacobian_opt, ode_opt);
end


convergenceStr = [];

% create unified odefun depending on mass usage
if NameValueArgs.isUseMassMatrix
    odefun_no_mass = @(tn, yn) [yn(dofNum+1:end); dynamicEquationNoMass(tn, yn(1:dofNum), yn(dofNum+1:end), Parameter)];
    odefun = odefun_no_mass;
else
    odefun_with_mass = @(tn, yn) [yn(dofNum+1:end); dynamicEquation(tn, yn(1:dofNum), yn(dofNum+1:end), Parameter)];
    odefun = odefun_with_mass;
end

% calculate with classic Runge-Kutta Method
if strcmp(NameValueArgs.calculateMethod, 'RK')
    q = zeros(dofNum, tNum); % for saving response
    dq = zeros(dofNum, tNum);
    equation = @(tn,yn,dyn) dynamicEquation(tn,yn,dyn,Parameter);
    % calculate response
    for iT = 1:tNum
        [q(:,iT), dq(:,iT)] = rungeKutta(equation, t(iT), yn, dyn, step);
        yn = q(:,iT);
        dyn = dq(:,iT);
        if any(isnan(dyn))
            convergenceStr = ['t=', num2str(t(iT)), 's non-convergent'];
            q = q(:,1:iT); dq = dq(:,1:iT); t = t(1:iT);
            break
        end
    end
elseif any(strcmp(NameValueArgs.calculateMethod, {'ode45','ode15s','ode23s'}))
    % select solver function handle
    switch NameValueArgs.calculateMethod
        case 'ode45', solver = @ode45;
        case 'ode15s', solver = @ode15s;
        case 'ode23s', solver = @ode23s;
    end
    % integrate
    [~, yode] = solver(odefun, t, [yn; dyn], ode_opt);
    yode = yode';
    q = yode(1:dofNum, :);
    dq = yode(dofNum+1:end, :);
else
    error('Error: wrong nameValue, please use RK or ode45 or ode15s or ode23s.')
end % end if NameValueArgs.calculateMethod

%%

%reduce data (re-sampling)
reduce_index = 1 : NameValueArgs.reduceInterval : size(q, 2); % index of the sampling data
if reduce_index(end) ~= size(q, 2)
   reduce_index = [reduce_index, size(q, 2)]; 
end
q  = q(:, reduce_index);% re_sampling
dq = dq(:,reduce_index);
t  = t(:,reduce_index);

%% 

% calculate the acceleration, rotational velocity and phase
if NameValueArgs.isPlotStatus
    % load constants
    Status   = Parameter.Status;
    shaftNum = Parameter.Shaft.amount;
    dofNum   = Parameter.Mesh.dofNum;
    vmax            = Status.vmax;
    duration        = Status.duration;
    acceleration    = Status.acceleration;
    isDeceleration  = Status.isDeceleration;
    vmin            = Status.vmin;
    ratio           = Status.ratio;
    ratio           = [1; ratio]; % the first shaft is basic

    
    % initial the status matrix
    omega 	= zeros(dofNum,length(t));
    domega  = omega;
    ddomega = omega;

    status = cell(shaftNum,1);

    if ~Status.isUseCustomize % using default status
        % calculate the status matrix
        for iShaft = 1:1:shaftNum
            iVmax           = vmax * ratio(iShaft);
            iVmin           = vmin * ratio(iShaft);
            iAcceleration   = acceleration*ratio(iShaft);
            status{iShaft} = rotationalStatus(t, iVmax, duration, iAcceleration,...
                                              isDeceleration, iVmin);
        end % end for iShaft
    else 
        % using customized status
        % initial
        for iShaft = 1:1:shaftNum
            status{iShaft} = zeros(3,length(t)); % initial
        end % end for iShaft
        % calculate status
        for iTime = 1:1:length(t)
            [ddomega_here, domega_here, omega_here] = Status.customize(t(iTime));
            for iShaft = 1:1:shaftNum
                status{iShaft}(1,iTime) = omega_here(iShaft);
                status{iShaft}(2,iTime) = domega_here(iShaft);
                status{iShaft}(3,iTime) = ddomega_here(iShaft);
            end % end for iShaft
        end % end for iTime
        
    end % end if


    % assemble the status matrix
    Node = Parameter.Mesh.Node;
    dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
    for iDof = 1:1:dofNum
        isBearing = Node( dofOnNodeNo(iDof) ).isBearing;
        if ~isBearing
            shaftNo = Node( dofOnNodeNo(iDof) ).onShaftNo;
            omega(iDof,:)   = status{shaftNo}(1,:);
            domega(iDof,:)  = status{shaftNo}(2,:);
            ddomega(iDof,:) = status{shaftNo}(3,:);
        end
    end


    % plot running status
    plotRunningStatus(t,status)
end



end

% sub function
function full_jacobian = jacobian_enable_mass_function(tn, yn, dyn, Parameter)

% get J21 and J22 without mass matrix
[J21, J22] = jacobianJ21J22(tn, yn, dyn, Parameter);

% construct full jacobian matrix
dofNum = Parameter.Mesh.dofNum;
full_jacobian = [sparse(dofNum, dofNum), speye(dofNum, dofNum);...
                 J21,                  J22];
end

% sub function
function full_jacobian = jacobian_disable_mass_function(tn, yn, dyn, Parameter)

% get J21 and J22 without mass matrix
[J21, J22] = jacobianJ21J22(tn, yn, dyn, Parameter);

% construct full jacobian matrix
dofNum = Parameter.Mesh.dofNum;
M = Parameter.Matrix.mass;
full_jacobian = [sparse(dofNum, dofNum), speye(dofNum, dofNum);...
                 M\J21,                    M\J22];

end