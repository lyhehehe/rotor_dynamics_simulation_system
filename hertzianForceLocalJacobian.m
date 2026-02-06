%% hertzianForceLocalJacobian - Calculate the local 2x2 tangent stiffness matrix for a bearing
%
% This function computes the partial derivatives of the Hertzian contact 
% force vector with respect to the relative radial displacements ($x$ and $y$). 
% It serves as the numerical kernel for assembling the global system Jacobian.
%
%% Syntax
%  K = hertzianForceLocalJacobian(x, y, thetai, thetao, nb, ri, ro, delta0, kHertz, n)
%
%% Description
% |hertzianForceLocalJacobian| implements the analytical derivative of the 
% nonlinear force-displacement relationship. For each rolling element in 
% contact, it calculates the contribution to the local stiffness matrix 
% based on the current cage phase and deformation state.
%
%% Input Arguments
% * |x|, |y| - Relative displacements between inner and outer races [m]
% * |thetai| - Inner race rotational phase [rad]
% * |thetao| - Outer race rotational phase [rad]
% * |nb| - Number of roller elements
% * |ri|, |ro| - Inner and outer race radii [m]
% * |delta0| - Radial clearance [m]
% * |kHertz| - Hertzian contact stiffness [N/m^n]
% * |n| - Contact force exponent (e.g., 1.5 for ball bearings)
%
%% Output Arguments
% * |K| - Local tangent stiffness matrix [2×2]:
%   $$ \mathbf{K} = \begin{bmatrix} k_{xx} & k_{xy} \\ k_{yx} & k_{yy} \end{bmatrix} $$
%
%% Mathematical Formulation
% The local stiffness matrix is the Jacobian of the force vector $\mathbf{f}$ 
% with respect to the displacement vector $\mathbf{q} = [x, y]^T$:
%
% 1. Local Deformation:
%    $\delta_k = x\cos\theta_k + y\sin\theta_k - \delta_0$
% 2. Derivative for a Single Roller (if $\delta_k > 0$):
%    $\frac{\partial f}{\partial \delta_k} = n \cdot \delta_k^{n-1}$
% 3. Stiffness Components:
%    $k_{xx} = k_{Hertz} \sum \left( n \delta_k^{n-1} \cos^2\theta_k \right)$
%    $k_{yy} = k_{Hertz} \sum \left( n \delta_k^{n-1} \sin^2\theta_k \right)$
%    $k_{xy} = k_{yx} = k_{Hertz} \sum \left( n \delta_k^{n-1} \cos\theta_k \sin\theta_k \right)$
%
%% Implementation Notes
% * **Symmetry**: The resulting matrix is symmetric ($k_{xy} = k_{yx}$), 
%   reflecting the conservative nature of the elastic contact force.
% * **Heaviside Condition**: Derivatives are only accumulated for rolling 
%   elements where the deformation $\delta_k$ is positive (in the contact zone).
%
%% See Also
% hertzianForceEq, hertzianForceJacobian, generateJacobian
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
function K = hertzianForceLocalJacobian(x, y, thetai, thetao, nb, ri, ro, delta0, kHertz, n)
%HERTZIANFORCELOCALJACOBIAN Calculates the local 2x2 stiffness matrix.

c1 = 2*pi/nb;
c2 = (thetao*ro + thetai*ri) / (ro + ri);

% Initialize stiffness accumulators
k11 = 0; % Sum for xx
k12 = 0; % Sum for xy (symmetric)
k22 = 0; % Sum for yy

for ik = 1:1:nb
    thetak = c1 * (ik-1) + c2; 
    c5 = cos(thetak);
    c6 = sin(thetak);
    deltak = x*c5 + y*c6 - delta0; 
    
    if deltak > 0
        % Derivative of deltak^n is: n * deltak^(n-1)
        stiff_val = n * (deltak^(n-1)); 
        
        % Reuse c5, c6 for matrix elements
        k11 = k11 + stiff_val * c5 * c5;
        k12 = k12 + stiff_val * c5 * c6;
        k22 = k22 + stiff_val * c6 * c6;
    end 
end

% Assemble local matrix
K = kHertz * [k11, k12; k12, k22];

end