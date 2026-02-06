%% generateJacobian - Automatic code generation for rotor system Jacobian matrices
%
% This function generates a MATLAB function file (jacobianJ21J22.m) that 
% computes the sub-matrices of the system Jacobian. These matrices are 
% essential for accelerating implicit numerical integration in rotor 
% dynamics simulations.
%
%% Syntax
%  generateJacobian(Parameter, calculateOmegaStr, loadMatrix1Str, loadMatrix2Str, processGNStr)
%
%% Description
% |generateJacobian| serves as a code-generation engine that ensures 
% mathematical consistency between the dynamic equations and their 
% derivatives. The generated function |jacobianJ21J22| computes:
% * The linear Jacobian components derived from system matrices (K, C, G, N)
% * The nonlinear Jacobian components from Hertzian contact models
% * The partial derivatives of user-defined custom forces
%
%% Input Arguments
% * |Parameter| - System configuration structure (used to check feature flags)
% * |calculateOmegaStr| - Code snippet for rotational kinematics calculation
% * |loadMatrix1Str| - Code snippet for loading primary system matrices (M, G, N, Q)
% * |loadMatrix2Str| - Code snippet for loading secondary matrices (K, C)
% * |processGNStr| - Code snippet for speed-dependent matrix scaling
%
%% Generated Function Details
% The generated |jacobianJ21J22.m| returns two sub-matrices, $J_{21}$ and $J_{22}$, 
% representing the partial derivatives of the acceleration-level equations:
%
% $$ J_{21} = \frac{\partial \ddot{q}}{\partial q}, \quad J_{22} = \frac{\partial \ddot{q}}{\partial \dot{q}} $$
%
% These are used to assemble the full state-space Jacobian matrix:
%
% $$ J_{full} = \begin{bmatrix} \mathbf{0} & \mathbf{I} \\ \mathbf{J}_{21} & \mathbf{J}_{22} \end{bmatrix} $$
%
%% Key Features
% 1. **Linear Consistency**: Automatically maps stiffness ($K$), damping ($C$), 
%    gyroscopic ($G$), and transient ($N$) matrices to their Jacobian counterparts.
% 2. **Hertzian Nonlinearity**: Incorporates |hertzianForceJacobian| to handle 
%    the state-dependent derivatives of rolling element contact.
% 3. **Custom Force Support**: Seamlessly integrates derivatives from 
%    |Parameter.Custom.jacobian| if provided by the user.
%
%% Implementation Notes
% * **String Injection**: This function relies on code snippets passed from 
%   |generateDynamicEquation| to ensure that kinematics and matrix loading 
%   logic are identical in both the EOM and Jacobian files.
% * **Performance**: The generated Jacobian utilizes sparse matrix operations 
%   to maintain efficiency in high-DOF systems.
%
%% Dependencies
% * |hertzianForceJacobian| - Derivative calculation for contact forces
% * Part of the code generation pipeline initiated by |generateDynamicEquation|
%
%% See Also
% generateDynamicEquation, calculateResponse, hertzianForceJacobian, establishModel
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
function generateJacobian(Parameter, calculateOmegaStr, loadMatrix1Str, loadMatrix2Str, processGNStr)

arguments
    Parameter 
    calculateOmegaStr 
    loadMatrix1Str 
    loadMatrix2Str
    processGNStr 
end
%% 

% check the exist of jacobian matrix and create .m
jacobianFileName = 'jacobianJ21J22.m';
if isfile(jacobianFileName)
    delete(jacobianFileName);
end

fJac = fopen(jacobianFileName,'w');


%% 

% write comments line
comments = [];


% write function start
functionStart = [...
"function [J21, J22] = jacobianJ21J22(tn, yn, dyn, Parameter)";...
" "...
];

fprintf(fJac,'%s\n',comments);
fprintf(fJac,'%s\n',functionStart);


%% 

% write status calculation
fprintf(fJac,'%s\n', calculateOmegaStr);


%%

% write matrix loading
fprintf(fJac,'%s\n', loadMatrix1Str);
fprintf(fJac,'%s\n', loadMatrix2Str);


%%

% write G N calculating
fprintf(fJac,'%s\n', processGNStr);

%% 

% write linear part output
linearOutput = [...
" ";...
"% calculate linear part of J21 and J22";...
"J21_lin = -(K - N);";...
"J22_lin = -(C - G);"
];
fprintf(fJac,'%s\n', linearOutput);


%%

% write Herzian nonlinear force part
% Hertzian contact force
if Parameter.ComponentSwitch.hasHertzianForce

    % write codes in dynamicEquation.m
    processHertzianForce = [
" ";...
"% calculate jacobian for Hertzian force";...
"dHerzian_dyn = hertzianForceJacobian(yn, omega, Parameter.Matrix.HerzianParameter, Parameter.Mesh.dofNum);";...
" ";...    
    ]; % write something
    fprintf(fJac,'%s\n', processHertzianForce);
    hertzJacobianJ21 = ' + dHerzian_dyn'; % write plus something
    hertzJacobianJ22 = '';% Herzian force does not contain dyn
else
    hertzJacobianJ21 = '';
    hertzJacobianJ22 = '';
end

%%

% write custom force part
if Parameter.ComponentSwitch.hasCustom
processCustomForce = [...
 " ";...
 "% calculate jacobian for customize force";...
 "[dfCustom_dyn, dfCustom_ddyn] = Parameter.Custom.jacobian(yn, dyn, tn, omega, domega, ddomega, Parameter);";...
 " ";...
];
    fprintf(fJac,'%s\n', processCustomForce);
    customJacobianJ21 = ' + dfCustom_dyn';
    customJacobianJ22 = ' + dfCustom_ddyn';
else
    customJacobianJ21 = '';
    customJacobianJ22 = '';
end

%% 

% write total Jacobian J21 and J22
totalJ21 = {...
 ' ';...
 '% total J21 ';...
 ['J21 = J21_lin', hertzJacobianJ21, customJacobianJ21, ';'];...
 ' ';...
 };
totalJ21 = cell2string(totalJ21);
fprintf(fJac,'%s\n', totalJ21);

totalJ22 = {...
 ' ';...
 '% total J22 ';...
 ['J22 = J22_lin', hertzJacobianJ22, customJacobianJ22, ';'];...
 ' ';...
 };
totalJ22 = cell2string(totalJ22);
fprintf(fJac,'%s\n', totalJ22);


%%

% write function end
functionEnd = [...
"end";...
" "...
];
fprintf(fJac,'%s\n',functionEnd);

%%

% close file

fclose(fJac);
end