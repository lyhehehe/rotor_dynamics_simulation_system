%% calculateUnbalanceResponse - Calculate steady-state unbalance response in the frequency domain
%
% This function solves the complex linear equation of motion at each operating
% condition to obtain the synchronous vibration amplitude and phase for every DOF.
% It supports multi-shaft systems with speed-dependent bearings and provides
% optional orbit trajectory plots.
%
%% Syntax
%  Response = calculateUnbalanceResponse(Parameter)
%  Response = calculateUnbalanceResponse(Parameter, speedMatrix)
%  Response = calculateUnbalanceResponse(Parameter, speedMatrix, Name=Value)
%
%% Description
% |calculateUnbalanceResponse| assembles the dynamic stiffness matrix
% H(omega) = -omega^2*M + i*omega*(C - G_total) + K at each operating condition
% and solves H*X = F for the complex displacement vector. The function:
% * Auto-generates a single-condition speed matrix from |Parameter.Status| if |speedMatrix| is omitted
% * Constructs the gyroscopic matrix dynamically for multi-shaft speed ratios
% * Updates speed-dependent bearing matrices via |updateSpeedDependentMatrices| when applicable
% * Superimposes unbalance forces from shafts sharing the same frequency before inversion
% * Separates the response by excitation source (one output slice per shaft)
% * Optionally plots orbit trajectories at specified nodes and conditions
%
%% Input Arguments
% * |Parameter| - System configuration structure containing:
%   * |Matrix|: Global matrices (mass, stiffness, damping, gyroscopic, unbalance)
%   * |Mesh|: Mesh data including |dofNum| and |dofInterval|
%   * |Status|: |vmax| (scalar, Shaft 1 max speed) and |ratio| for multi-shaft speed scaling
%   * |Shaft|: Shaft properties including |amount|
% * |speedMatrix| - Operating speeds [rad/s]; each column is one condition [shaftNum×numSpeeds double]
%   Pass |[]| to use a single default condition from |Parameter.Status.vmax| and |ratio|
%
%% Name-Value Arguments
% * |isPlot|           - Plot orbit trajectories (default: false) [logical]
% * |plotNodeID|       - Node IDs to plot (default: 1) [1×P double]
% * |plotConditionIdx| - Speed condition indices to plot (default: 1) [1×Q double]
%
%% Output Arguments
% * |Response| - Complex steady-state displacement [m] [dofNum×numSpeeds×shaftNum double]
%   Response(:, i, k) is the vibration at condition i excited exclusively by
%   unbalance forces from Shaft k (at Shaft k's rotation frequency)
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function Response = calculateUnbalanceResponse(Parameter, speedMatrix, options)

    arguments
        Parameter struct
        speedMatrix double = []
        options.isPlot (1,1) logical = false
        options.plotNodeID (1,:) double = 1
        options.plotConditionIdx (1,:) double = 1
    end

    % 1. Speed Matrix Generation & Validation
    shaftNum = Parameter.Shaft.amount;
    
    % Auto-generate speed matrix if not provided (or passed as [])
    if isempty(speedMatrix)
        vmax = Parameter.Status.vmax; 
        speedMatrix = zeros(shaftNum, 1);
        speedMatrix(1, 1) = vmax; % Reference shaft (Shaft 1)
        
        if shaftNum > 1
            speedRatio = Parameter.Status.ratio; % Length should be shaftNum - 1
            for iShaft = 2:shaftNum
                speedMatrix(iShaft, 1) = vmax * speedRatio(iShaft-1);
            end
        end
    end
    
    if size(speedMatrix, 1) ~= shaftNum
        error('The number of rows in speedMatrix must equal the number of shafts.');
    end
    numSpeeds = size(speedMatrix, 2);

    % 2. Extract System Matrices and Mesh Info
    M = Parameter.Matrix.mass;
    C = Parameter.Matrix.damping;
    K = Parameter.Matrix.stiffness;
    G_base = Parameter.Matrix.gyroscopic;
    
    dofNum = Parameter.Mesh.dofNum;
    dofInterval = Parameter.Mesh.dofInterval;
    Node = Parameter.Mesh.Node;
    unbalanceData = Parameter.Matrix.unbalance; % [NodeID, ShaftID, Magnitude, Phase]
    
    % Get shaft DOF range for gyroscopic matrix assembly
    shaftDof = zeros(shaftNum, 2);
    for iShaft = 1:shaftNum
        IShaftNode = Node([Node.onShaftNo] == iShaft & [Node.isBearing] == false);
        startNode = min([IShaftNode.name]);
        endNode = max([IShaftNode.name]);
        shaftDof(iShaft,:) = [dofInterval(startNode,1), dofInterval(endNode,2)];
    end

    % 3. Pre-assemble Unbalance Base Vector (U_base)
    U_base = zeros(dofNum, shaftNum);
    for i = 1:size(unbalanceData, 1)
        nodeID = unbalanceData(i, 1);
        shaftID = unbalanceData(i, 2);
        mag = unbalanceData(i, 3);
        phase = unbalanceData(i, 4); % Phase is strictly in radians
        
        nodeStartDof = dofInterval(nodeID, 1);
        dofX = nodeStartDof;
        dofY = nodeStartDof + 1;
        
        % Right-hand rule convention: Y lags X by 90 degrees
        U_base(dofX, shaftID) = U_base(dofX, shaftID) + mag * exp(1i * phase);
        U_base(dofY, shaftID) = U_base(dofY, shaftID) - 1i * mag * exp(1i * phase);
    end

    % 4. Pre-allocate Response Output Matrix
    Response = zeros(dofNum, numSpeeds, shaftNum);
    speedTol = 1e-5; 

    % 5. Main Calculation Loop
    for iSpeed = 1:numSpeeds
        currentSpeeds = speedMatrix(:, iSpeed);
        
        % Assemble Total Gyroscopic Matrix
        G_total = G_base; 
        for iShaft = 1:shaftNum
            rng = shaftDof(iShaft,1) : shaftDof(iShaft,2);
            G_total(rng, rng) = currentSpeeds(iShaft) * G_total(rng, rng);
        end
        
        % Update Speed-Dependent Matrices ONCE per steady-state
        if isfield(Parameter, 'ComponentSwitch') && ...
           isfield(Parameter.ComponentSwitch, 'hasSpeedDependentBearing') && ...
           Parameter.ComponentSwitch.hasSpeedDependentBearing
            
            % Pass the full 'currentSpeeds' vector so each bearing knows its shaft's speed
            [K_active, C_active] = updateSpeedDependentMatrices(...
                Parameter.SpeedDependentBearing, K, C, currentSpeeds, dofInterval);
        else
            K_active = K;
            C_active = C;
        end
        
        activeShafts = find(currentSpeeds > speedTol);
        if isempty(activeShafts), continue; end
        
        uniqueSpeeds = unique(currentSpeeds(activeShafts));
        
        % Excitation Frequency Loop (omega)
        for j = 1:length(uniqueSpeeds)
            omega = uniqueSpeeds(j);
            sameSpeedShafts = find(abs(currentSpeeds - omega) < speedTol);
            
            U_combined = sum(U_base(:, sameSpeedShafts), 2);
            F_total = U_combined * (omega^2);
            
            if max(abs(F_total)) > 1e-12
                % Use the constant K_active and C_active for this steady-state condition
                H = -(omega^2) * M + 1i * omega * (C_active - G_total) + K_active;
                X = H \ F_total;
                
                for k = 1:length(sameSpeedShafts)
                    shaftIdx = sameSpeedShafts(k);
                    Response(:, iSpeed, shaftIdx) = X;
                end
            end
        end
    end
    
    % 6. Post-Processing: Orbit Plotting Integration
    if options.isPlot
        for pNode = options.plotNodeID
            if pNode > size(dofInterval, 1) || pNode < 1
                warning('Node %d does not exist. Skipping plot.', pNode);
                continue;
            end
            for pCond = options.plotConditionIdx
                if pCond > numSpeeds || pCond < 1
                    warning('Condition index %d is out of bounds. Skipping plot.', pCond);
                    continue;
                end
                plotOrbitLocal(Response, Parameter, speedMatrix, pNode, pCond);
            end
        end
    end
end

%% Local Function for Orbit Plotting
function plotOrbitLocal(Response, Parameter, speedMatrix, targetNodeID, speedIndex)
% Extracts physical time-domain responses and plots the total orbit trajectory
    shaftNum = Parameter.Shaft.amount;
    dofInterval = Parameter.Mesh.dofInterval;
    
    dofX = dofInterval(targetNodeID, 1);
    dofY = dofX + 1;
    
    currentSpeeds = speedMatrix(:, speedIndex);
    activeSpeeds = currentSpeeds(currentSpeeds > 1e-5);
    
    if isempty(activeSpeeds)
        warning('Node %d Condition %d: All shafts are at zero speed. No orbit to plot.', targetNodeID, speedIndex);
        return;
    end
    
    % Adaptive time window: enough cycles of the slowest shaft
    minSpeed = min(activeSpeeds);
    maxTime = 10 * (2 * pi / minSpeed); 
    
    % Adaptive sample rate: high enough resolution for the fastest shaft
    maxSpeed = max(activeSpeeds);
    dt = (2 * pi / maxSpeed) / 150; % 150 points per fastest revolution
    t = 0:dt:maxTime;
    
    x_total = zeros(size(t));
    y_total = zeros(size(t));
    
    % Superimpose physical (real) time-domain displacements from all excitation sources
    for k = 1:shaftNum
        omega = currentSpeeds(k);
        if omega > 1e-5
            X_complex = Response(dofX, speedIndex, k);
            Y_complex = Response(dofY, speedIndex, k);
            
            x_total = x_total + real(X_complex * exp(1i * omega * t));
            y_total = y_total + real(Y_complex * exp(1i * omega * t));
        end
    end
    
    % Plot Generation
    figure('Name', sprintf('Orbit of Node %d', targetNodeID));
    plot(x_total, y_total, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2); % Standard MATLAB blue
    hold on;
    
    % Mark the starting coordinate
    plot(x_total(1), y_total(1), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Start (t=0)');
    
    grid on;
    xlabel('X Displacement (m)', 'FontSize', 11);
    ylabel('Y Displacement (m)', 'FontSize', 11);
    
    % Maintained your exact bounding logic
    all_data = [x_total; y_total];
    max_val = max(abs(all_data), [], 'all'); 
    limit_val = max_val * 1.1; 
    
    xlim([-limit_val, limit_val]);
    ylim([-limit_val, limit_val]);
    axis square;
    
    % Convert speeds to RPM for the title string
    rpmStr = sprintf('%.0f ', currentSpeeds * 30/pi);
    title({sprintf('Orbit Trajectory of Node %d', targetNodeID), ...
           sprintf('Condition %d | Speeds: [%s] RPM', speedIndex, strtrim(rpmStr))}, ...
          'FontSize', 12, 'FontWeight', 'bold');
          
    xline(0, 'k--', 'HandleVisibility', 'off');
    yline(0, 'k--', 'HandleVisibility', 'off');
    hold off;
end