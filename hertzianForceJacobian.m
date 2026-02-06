%% hertzianForceJacobian - Calculate the tangent stiffness matrix of Hertzian contact forces
%
% This function computes the Jacobian matrix (partial derivative of Hertzian 
% forces with respect to displacement, $\partial F_{hertz} / \partial q$) 
% for rotor-bearing systems. It is primarily used by implicit ODE solvers 
% to improve convergence and stability.
%
%% Syntax
%  JHertz = hertzianForceJacobian(qn, omega, HerzianParameter, dofNum)
%
%% Description
% |hertzianForceJacobian| constructs a global sparse tangent stiffness matrix 
% by assembling local Jacobian contributions from each Hertzian bearing. 
% The function:
% * Computes local 2x2 Jacobian kernels for each contact pair
% * Assembles 4x4 interaction blocks for inner/outer race coupling
% * Handles ground-connected bearings via index filtering
% * Utilizes efficient sparse triplet format (row, col, val) for assembly
%
%% Input Arguments
% * |qn| - Current system displacement vector [dofNum × 1]
% * |omega| - Angular phase vector for each shaft [shaftNum × 1]
% * |HerzianParameter| - Consolidated contact parameters containing:
%   * |nb|, |ri|, |ro|, |delta0|, |kHertz|, |n|: Bearing geometric/material constants
%   * |omegaiNo|, |omegaoNo|: Shaft speed indices
%   * |hertzDof|: DOF mapping indices [N × 2]
%   * |hertzianNum|: Number of active Hertzian contact elements
% * |dofNum| - Total number of system degrees of freedom
%
%% Output Arguments
% * |JHertz| - Global sparse Jacobian matrix [dofNum × dofNum]
%
%% Mathematical Formulation
% For a bearing connecting an inner race ($i$) and an outer race ($o$), the 
% force depends on relative displacement $\Delta q = q_i - q_o$. The 
% contribution to the global Jacobian is structured as a 4x4 block:
%
% $$ \mathbf{K}_{block} = \begin{bmatrix} 
% \frac{\partial \mathbf{F}_i}{\partial \mathbf{q}_i} & \frac{\partial \mathbf{F}_i}{\partial \mathbf{q}_o} \\
% \frac{\partial \mathbf{F}_o}{\partial \mathbf{q}_i} & \frac{\partial \mathbf{F}_o}{\partial \mathbf{q}_o}
% \end{bmatrix} = \begin{bmatrix} 
% -\mathbf{K}_{loc} & \mathbf{K}_{loc} \\
% \mathbf{K}_{loc} & -\mathbf{K}_{loc}
% \end{bmatrix} $$
%
% where $\mathbf{K}_{loc}$ is the 2x2 local derivative $\frac{\partial \mathbf{F}}{\partial \Delta \mathbf{q}}$.
%
%% Implementation Notes
% * **Sparse Assembly**: To maximize performance, the matrix is built using 
%   pre-allocated triplet vectors which are passed to the |sparse| constructor.
% * **Ground Adaptation**: Bearings connected to the ground have outer 
%   race indices pointing to virtual DOFs ($> dofNum$), which are 
%   filtered out during the final matrix construction.
%
%% Dependencies
% * |hertzianForceLocalJacobian| - Local 2x2 derivative kernel
%
%% See Also
% hertzianForce, hertzianForceLocalJacobian, generateJacobian, calculateResponse
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
function JHertz = hertzianForceJacobian(qn, omega, HerzianParameter, dofNum)
%HERTZIANFORCEJACOBIAN Calculates the Jacobian matrix (tangent stiffness) of Hertzian contact forces.

% Load parameters
nb = HerzianParameter.nb;
ri = HerzianParameter.ri;
ro = HerzianParameter.ro;
delta0 = HerzianParameter.delta0;
kHertz = HerzianParameter.kHertz;
n = HerzianParameter.n;
omegaiNo = HerzianParameter.omegaiNo;
omegaoNo = HerzianParameter.omegaoNo;
hertzDof = HerzianParameter.hertzDof;
hertzianNum = HerzianParameter.hertzianNum;

% Pre-processing input vectors for ground adaptation
omega = [omega, 0]; 
qn = [qn; 0; 0]; 

% Pre-allocate arrays for sparse matrix triplets
% Each bearing contributes to a 4x4 block (16 elements max)
maxNz = hertzianNum * 16;
rows = zeros(maxNz, 1);
cols = zeros(maxNz, 1);
vals = zeros(maxNz, 1);
count = 1;

for iHertz = 1:1:hertzianNum
    % Extract indices
    idxInX = hertzDof(iHertz,1);
    idxInY = hertzDof(iHertz,1)+1;
    idxOutX = hertzDof(iHertz,2);
    idxOutY = hertzDof(iHertz,2)+1;
    
    % Calculate relative displacement (Inner - Outer)
    x = qn(idxInX) - qn(idxOutX);
    y = qn(idxInY) - qn(idxOutY);
    
    % Call local Jacobian kernel
    K_loc = hertzianForceLocalJacobian(x, y, omega(omegaiNo(iHertz)), omega(omegaoNo(iHertz)), ...
        nb(iHertz), ri(iHertz), ro(iHertz), delta0(iHertz), kHertz(iHertz), n(iHertz));
    
    % Assemble 4x4 block into triplets
    % The structure is: [ -K_loc,  K_loc; 
    %                      K_loc, -K_loc ]
    
    % Define Global Indices for this block
    gIdx = [idxInX, idxInY, idxOutX, idxOutY];
    
    % Iterate through the 2x2 local matrix to fill the 4x4 global block
    for r = 1:2
        for c = 1:2
            val = K_loc(r,c);   
            % Block (1,1): Inner-Inner (Partial Force_In / Partial Disp_In) -> -K_loc
            rows(count) = gIdx(r); cols(count) = gIdx(c); vals(count) = -val; count = count + 1;
            
            % Block (1,2): Inner-Outer (Partial Force_In / Partial Disp_Out) -> +K_loc
            rows(count) = gIdx(r); cols(count) = gIdx(c+2); vals(count) = val; count = count + 1;
            
            % Block (2,1): Outer-Inner (Partial Force_Out / Partial Disp_In) -> +K_loc
            rows(count) = gIdx(r+2); cols(count) = gIdx(c); vals(count) = val; count = count + 1;
            
            % Block (2,2): Outer-Outer (Partial Force_Out / Partial Disp_Out) -> -K_loc
            rows(count) = gIdx(r+2); cols(count) = gIdx(c+2); vals(count) = -val; count = count + 1;
        end
    end
end

% Filter out valid entries (remove unused pre-allocated slots)
validIdx = 1:(count-1);
rows = rows(validIdx);
cols = cols(validIdx);
vals = vals(validIdx);

% Filter out Ground DOFs (indices > dofNum)
mask = (rows <= dofNum) & (cols <= dofNum);
JHertz = sparse(rows(mask), cols(mask), vals(mask), dofNum, dofNum);

end