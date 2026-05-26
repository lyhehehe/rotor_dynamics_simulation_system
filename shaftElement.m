%% shaftElement - Generate finite element matrices for Timoshenko beam shaft elements
%
% This function computes mass, stiffness, gyroscopic, and transient matrices 
% along with gravity forces and unbalance data for Timoshenko beam elements.
% It supports dual-geometry definitions (physical vs. stiffness diameters) 
% to model stepped shafts or stiffness-equivalent diameters accurately.
%
%% Syntax
%   [Me, Ke, Ge, Ne, Fge, UnbalInfo] = shaftElement(ElemProp)
%
%% Description
% |shaftElement| calculates finite element matrices for shaft segments using 
% Timoshenko beam theory, accounting for:
% * Shear deformation effects (using stiffness geometry)
% * Rotary inertia (using physical geometry)
% * Gyroscopic moments
% * Transient vibration
% * Gravity loading
% * Distributed mass unbalance
%
%% Input Arguments
% * |ElemProp| - Element properties structure with fields:
%   * |dofOfEachNodes|: DOF per node [scalar]
%   * |Length|: Element length [m]
%   * |density|: Material density [kg/m³]
%   * |elasticModulus|: Elastic modulus [Pa]
%   * |poissonRatio|: Poisson's ratio
%   * |outerRadius|: Physical outer radius [m] (for Mass/Inertia)
%   * |innerRadius|: Physical inner radius [m] (for Mass/Inertia)
%   * |outerRadiusStiff|: Effective outer radius [m] (for Stiffness)
%   * |innerRadiusStiff|: Effective inner radius [m] (for Stiffness)
%   * |eccentricity|: Mass eccentricity [m]
%   * |eccentricityPhase|: Phase of eccentricity [rad]
%
%% Output Arguments
% * |Me| - Mass matrix [8×8] (based on physical geometry)
% * |Ke| - Stiffness matrix [8×8] (based on stiffness geometry)
% * |Ge| - Gyroscopic matrix [8×8]
% * |Ne| - Transient matrix [8×8]
% * |Fge| - Gravity force vector [8×1]
% * |UnbalInfo| - Unbalance structure:
%       .Magnitude - Static unbalance moment [kg·m] ($m_{elem} \cdot e$)
%       .Phase     - Phase angle [rad]
%
%% Physical Parameters
% * Geometry (Mass/Gyro):
%   $A_{geo} = \pi(R_{geo}^2 - r_{geo}^2)$, $I_{geo} = \frac{\pi}{4}(R_{geo}^4 - r_{geo}^4)$
% * Stiffness (Deformation):
%   $A_{stiff} = \pi(R_{stiff}^2 - r_{stiff}^2)$, $I_{stiff} = \frac{\pi}{4}(R_{stiff}^4 - r_{stiff}^4)$
% * Shear Deformation ($\phi_s$):
%   Calculated using $I_{stiff}$ and $A_{shear}$ derived from stiffness radii.
%
%% Matrix Formulation
% 1. Mass Matrix (Me):
%    * Uses $\rho_L = \rho \cdot A_{geo}$ and $I_{geo}$
%    * Includes Translational ($M_T$) and Rotational ($M_R$) components
% 2. Stiffness Matrix (Ke):
%    * Uses $E \cdot I_{stiff}$ and shear coefficient derived from stiffness radii
% 3. Gyroscopic Matrix (Ge):
%    * Derived from $N_e$ using polar inertia ($2 \cdot I_{geo}$)
% 4. Gravity Vector (Fge):
%    * $F_{ge} = -(\rho \cdot A_{geo} \cdot l) \cdot g / 2$ at translational DOFs
% 5. Unbalance:
%    * $U = (\rho \cdot A_{geo} \cdot l) \cdot e$
%
%% Implementation Notes
% 1. Element DOF Ordering:
%    [x1, y1, θx1, θy1, x2, y2, θx2, θy2]^T
% 2. Dual-Radius Model:
%    Separates physical dimensions (mass) from stiffness dimensions. 
%    Useful for simplifying complex geometries into beams while retaining 
%    correct inertial properties.
% 3. Hollow Shaft Handling:
%    Supports annular cross-sections for both physical and stiffness definitions.
%
%% Example
% % Define element properties
% elemProps = struct(...
%     'Length', 0.1, ...
%     'outerRadius', 0.05, 'innerRadius', 0.03, ...       % Physical
%     'outerRadiusStiff', 0.05, 'innerRadiusStiff', 0.03, ... % Stiffness
%     'density', 7850, 'elasticModulus', 210e9, ...
%     'poissonRatio', 0.3, ...
%     'eccentricity', 1e-4, 'eccentricityPhase', 0);
% 
% % Generate element matrices
% [Me, Ke, Ge, Ne, Fge, Unbal] = shaftElement(elemProps);
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


function [Me, Ke, Ge, Ne, Fge, UnbalInfo] = shaftElement(ElemProp)

%% 1. Extract and Calculate Properties

% A. Geometric Properties (For Mass M, Gyroscopic G, Gravity)
r_geo = ElemProp.innerRadius;
R_geo = ElemProp.outerRadius;
A_geo = pi*R_geo^2 - pi*r_geo^2;        % Geometric Area
I_geo = pi/4 * (R_geo^4 - r_geo^4);     % Geometric Area Moment
% Ip_geo = 2 * I_geo;                     % Polar Moment of Inertia

% B. Stiffness Properties (For Stiffness K)
r_stiff = ElemProp.innerRadiusStiff;
R_stiff = ElemProp.outerRadiusStiff;
% Note: Shear area (As) is typically related to the cross-section providing stiffness.
% Hence we use stiffness radii for As and I.
A_stiff = pi*R_stiff^2 - pi*r_stiff^2;
I_stiff = pi/4 * (R_stiff^4 - r_stiff^4);

% C. Material & Dimensions
l   = ElemProp.Length;
E   = ElemProp.elasticModulus;
mu  = ElemProp.poissonRatio;
rho = ElemProp.density;

% D. Timoshenko Beam Coefficients
% Shear coefficient calculation (using stiffness geometry)
As1 = (7 + 6*mu) / (6 * (1 + mu)); 
As2 = (20 + 12*mu) / (7 + 6*mu); 
if (R_stiff^2 + r_stiff^2) == 0
    As3 = 0; % Avoid div by zero for phantom elements
else
    As3 = ((R_stiff * r_stiff) / (R_stiff^2 + r_stiff^2))^2;
end
As = A_stiff / (As1 * (1 + As2 * As3)); % Shear Area

% Shear deformation parameter (Phi)
% Using Stiffness I and Stiffness Area
% phis = 12 * E * I_stiff / (l^2 * (As * E / (2*(1+mu)))); 
% Simplified: phis = 24*I*(1+mu) / (As*l^2)
phis = 24 * I_stiff * (1 + mu) / (As * l^2);

% E. Linear Density (Physical Mass)
% Mass depends on geometric volume and density
rhoL = rho * A_geo; 

%% 2. Mass Matrix (Translation) - Uses Physical Mass (rhoL)
coefficient = rhoL * l / (1+phis)^2;
MT1 = 13/35 + (7/10)*phis + (1/3)*phis^2;
MT2 = l^2 * ( 1/105 +(1/60)*phis + (1/120)*phis^2 );
MT3 = 9/70 + (3/10)*phis + (1/6)*phis^2;
MT4 = l * ( 11/210 + (11/120)*phis + (1/24)*phis^2 );
MT5 = l * ( 13/420 + (3/40)*phis + (1/24)*phis^2 );
MT6 = (-1)*l^2 * ( 1/140 + (1/60)*phis + (1/120)*phis^2 );

MT = [ MT1,    0,    0,    0,    0,    0,    0,    0;...
         0,  MT1,    0,    0,    0,    0,    0,    0;...
         0, -MT4,  MT2,    0,    0,    0,    0,    0;...
       MT4,    0,    0,  MT2,    0,    0,    0,    0;...
       MT3,    0,    0,  MT5,  MT1,    0,    0,    0;...
         0,  MT3, -MT5,    0,    0,  MT1,    0,    0;...
         0,  MT5,  MT6,    0,    0,  MT4,  MT2,    0;...
      -MT5,    0,    0,  MT6, -MT4,    0,    0,  MT2];
MT = coefficient * triangular2symmetric(MT);

%% 3. Mass Matrix (Rotation) - Uses Physical Inertia
% Note: Using I_geo for rotational inertia
coefficient = rhoL * I_geo / ( l * (1+phis)^2 * A_geo ); 

MR1 = 6/5;
MR2 = l^2 * ( 2/15 + (1/6)*phis + (1/3)*phis^2);
MR3 = l^2 * ( -1/30 - (1/6)*phis + (1/6)*phis^2 );
MR4 = l * (1/10 - (1/2)*phis);

MR = [ MR1,    0,    0,    0,    0,    0,    0,    0;...
         0,  MR1,    0,    0,    0,    0,    0,    0;...
         0, -MR4,  MR2,    0,    0,    0,    0,    0;...
       MR4,    0,    0,  MR2,    0,    0,    0,    0;...
      -MR1,    0,    0, -MR4,  MR1,    0,    0,    0;...
         0, -MR1,  MR4,    0,    0,  MR1,    0,    0;...
         0, -MR4,  MR3,    0,    0,  MR4,  MR2,    0;...
       MR4,    0,    0,  MR3, -MR4,    0,    0,  MR2];
MR = coefficient * triangular2symmetric(MR);

% Total Element Mass Matrix
Me = MT + MR;

%% 4. Circulatory Matrix (N) - Uses Physical Inertia
% Related to rotary inertia effect
coefficient = rhoL * I_geo / ( 15 * l * (1+phis)^2 * A_geo );

N1 = 36;
N2 = 3*l - 15*l*phis;
N3 = l^2 + 5*l^2*phis - 5*l^2*phis^2;
N4 = 4*l^2 + 5*l^2*phis + 10*l^2*phis^2;

Ne = [  0, -N1,  N2,   0,   0,  N1,  N2,   0;...
        0,   0,   0,   0,   0,   0,   0,   0;...
        0,   0,   0,   0,   0,   0,   0,   0;...
        0, -N2,  N4,   0,   0,  N2, -N3,   0;...
        0   N1, -N2,   0,   0, -N1, -N2,   0;...
        0,   0,   0,   0,   0,   0,   0,   0;...
        0,   0,   0,   0,   0,   0,   0,   0;...
        0, -N2, -N3,   0,   0,  N2,  N4,   0 ];
Ne = coefficient * Ne;

%% 5. Gyroscopic Matrix (G)
Ge = Ne - Ne';

%% 6. Stiffness Matrix (K) - Uses Stiffness Properties (I_stiff)
coefficient = E * I_stiff / ( l^3 * (1+phis) );

K1 = 12;
K2 = l^2 * ( 4 + phis );
K3 = l^2 * ( 2 - phis );
K4 = 6*l;

Ke = [ K1,   0,   0,   0,   0,   0,   0,   0;...
        0,  K1,   0,   0,   0,   0,   0,   0;...
        0, -K4,  K2,   0,   0,   0,   0,   0;...
       K4,   0,   0,  K2,   0,   0,   0,   0;...
      -K1,   0,   0, -K4,  K1,   0,   0,   0;...
        0, -K1,  K4,   0,   0,  K1,   0,   0;...
        0, -K4,  K3,   0,   0,  K4,  K2,   0;...
       K4,   0,   0,  K3, -K4,   0,   0,  K2 ];
Ke = coefficient * triangular2symmetric(Ke);

%% 7. Gravity Vector
% Based on physical mass
m_elem = rhoL * l; 
FgeTotal = m_elem * 9.8; 
Fge = [0; -FgeTotal/2; 0; 0; 0; -FgeTotal/2; 0; 0];

%% 8. Unbalance Information
% Calculate total unbalance moment for this element
% U = mass * eccentricity
U_mag = m_elem * ElemProp.eccentricity;
U_phase = ElemProp.eccentricityPhase;

UnbalInfo.Magnitude = U_mag;
UnbalInfo.Phase = U_phase;

end