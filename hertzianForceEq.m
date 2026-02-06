%% hertzianForceEq - Calculate Hertzian contact forces in roller bearings
%
% This function computes the nonlinear Hertzian contact forces between 
% rollers and races in rolling element bearings based on relative 
% displacements and angular phases of the races.
%
%% Syntax
%  f = hertzianForceEq(x, y, thetai, thetao, nb, ri, ro, delta0, kHertz, n)
%
%% Description
% |hertzianForceEq| calculates the dynamic contact forces in rolling 
% element bearings using Hertzian contact theory. The function:
% * Models roller-race interactions
% * Accounts for cage angular position (phase)
% * Handles multiple roller elements
% * Implements nonlinear deformation effects
%
%% Input Arguments
% * |x|, |y| - Displacements in x/y directions [m]
% * |thetai| - Inner race rotational phase [rad]
% * |thetao| - Outer race rotational phase [rad]
% * |nb| - Number of roller elements [integer]
% * |ri| - Inner race radius [m]
% * |ro| - Outer race radius [m]
% * |delta0| - Initial roller-race clearance [m]
% * |kHertz| - Hertzian contact stiffness [N/m^n]
% * |n| - Contact force exponent [dimensionless]
%   * Typical values: 1.5 (ball bearings), 10/9 (roller bearings)
%
%% Output Arguments
% * |f| - Contact force vector [2×1]:
%   * [fx; fy] in [N]
%
%% Physical Formulation
% 1. Cage Phase (Angular Position):
%    $\theta_c = \frac{\theta_o r_o + \theta_i r_i}{r_o + r_i}$
% 2. Roller Position:
%    $\theta_k = \frac{2\pi}{n_b}(k-1) + \theta_c$
% 3. Local Deformation:
%    $\delta_k = x\cos\theta_k + y\sin\theta_k - \delta_0$
% 4. Contact Force:
%    $f_x = k_{Hertz} \sum\limits_{k=1}^{n_b} \delta_k^n \cos\theta_k \cdot H(\delta_k)$
%    $f_y = k_{Hertz} \sum\limits_{k=1}^{n_b} \delta_k^n \sin\theta_k \cdot H(\delta_k)$
%    where $H(\delta_k) = \begin{cases} 
%        1 & \delta_k > 0 \\
%        0 & \text{otherwise}
%    \end{cases}$
%
%% Implementation Details
% 1. Cage Motion:
%   * Computes the average orbital phase from race angular positions
% 2. Roller Loop:
%   * Iterates through all roller elements
%   * Calculates angular position $\theta_k$ for each roller
% 3. Deformation Check:
%   * Only considers positive deformations (contact condition)
% 4. Force Accumulation:
%   * Sums vector components with nonlinear stiffness
%
%% Application Notes
% * Bearings Supported:
%   * Ball bearings (n=3/2)
%   * Cylindrical roller bearings (n=10/9)
%   * Tapered roller bearings (n=10/9)
% * Force Direction:
%   * Positive forces indicate roller pushing race outward
% * Deformation Limitations:
%   * Only elastic deformations within Hertzian theory limits
%
%% Example
% % Calculate contact force in ball bearing
% displacement = [0.1e-3; -0.05e-3]; % x/y displacements [m]
% phi_i = pi/2;             % Inner race phase [rad]
% phi_o = 0;                % Outer race phase (ground) [rad]
% balls = 8;                % Number of balls
% ri = 0.015; ro = 0.025;   % Race radii [m]
% clearance = 1e-4;         % Initial clearance [m]
% stiff = 1.5e8;            % Contact stiffness [N/m^1.5]
% exp = 1.5;                % Ball bearing exponent
%
% % Calculate force
% f = hertzianForceEq(displacement(1), displacement(2), ...
%                     phi_i, phi_o, balls, ri, ro, clearance, stiff, exp);
%
%% See Also
% hertzianForce, hertzianForceJacobian, collectHerzianParameter
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function f = hertzianForceEq(x, y, thetai, thetao, nb, ri, ro, delta0, kHertz, n)

c1 = 2*pi/nb;
c2 = (thetao*ro + thetai*ri) / (ro + ri); % rotational angle of the cage;
c3 = 0;
c4 = 0;
for ik = 1:1:nb
    thetak = c1 * (ik-1) + c2; % the angular position of the k-th roller in the bearing
    c5 = cos(thetak);
    c6 = sin(thetak);
    deltak = x*c5 + y*c6 - delta0; % the relative contact deformation between the k-th roller and the inner race of the bearing
    if deltak>0
        c7 = deltak^n;
        c3 = c3 + c7 * c5;
        c4 = c4 + c7 * c6;
    end % end if
end
f = kHertz * [c3; c4];

end