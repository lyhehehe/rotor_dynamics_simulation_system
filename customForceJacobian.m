%% customForceJacobian - Calculate Jacobian matrices of custom external forces
%
% This function provides the partial derivatives of user-defined external 
% forces with respect to the system states (displacement and velocity), 
% which is essential for implicit ODE solvers (e.g., ode15s, ode23s).
%
%% Syntax
%  [dfCustomd_dyn, dfCustomd_ddyn] = customForceJacobian(yn, dyn, tn, omega, domega, ddomega, Parameter)
%
%% Description
% |customForceJacobian| computes the Jacobian matrices for custom force 
% components. This is a template function that should be populated with 
% the analytical or numerical derivatives of the forces defined in the 
% system. The function returns:
% * The partial derivative of forces with respect to displacement (Stiffness-like)
% * The partial derivative of forces with respect to velocity (Damping-like)
%
%% Input Arguments
% * |yn| - Displacement vector at current time step [dofNum × 1]
% * |dyn| - Velocity vector at current time step [dofNum × 1]
% * |tn| - Current simulation time [s]
% * |omega| - Angular phase vector for each shaft [shaftNum × 1]
% * |domega| - Angular velocity vector for each shaft [shaftNum × 1]
% * |ddomega| - Angular acceleration vector for each shaft [shaftNum × 1]
% * |Parameter| - System configuration structure containing Mesh and Shaft data
%
%% Output Arguments
% * |dfCustomd_dyn| - Jacobian matrix with respect to displacement [dofNum × dofNum]
%   (Partial derivative: ∂F_custom / ∂q)
% * |dfCustomd_ddyn| - Jacobian matrix with respect to velocity [dofNum × dofNum]
%   (Partial derivative: ∂F_custom / ∂dq)
%
%% Implementation Note
% To ensure performance, the output matrices should be returned as sparse 
% matrices. These matrices are used by the solver to construct the full 
% system Jacobian:
% J_full = [ 0,  I;  M \ (J21),  M \ (J22) ]
% where J21 and J22 incorporate these custom force derivatives.
%
%% See Also
% calculateResponse, dynamicEquation, jacobianJ21J22
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
% 
function [dfCustomd_dyn, dfCustomd_ddyn] = customForceJacobian(yn, dyn, tn, omega, domega, ddomega, Parameter)

dof_num = Parameter.Mesh.dofNum;

dfCustomd_dyn = sparse(dof_num, dof_num);
dfCustomd_ddyn = dfCustomd_dyn;

end