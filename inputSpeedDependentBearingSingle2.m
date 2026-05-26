%% inputSpeedDependentBearingSingle2 - Configure speed-dependent bearing parameters
%
% This function configures bearings whose stiffness and damping matrices are
% tabulated as functions of shaft speed. It represents linearized bearing
% coefficients (e.g., derived from Hertzian contact analysis) for use with
% the frequency-domain unbalance response solver.
%
%% Syntax
%  Parameter = inputSpeedDependentBearingSingle2(InitialParameter)
%
%% Description
% |inputSpeedDependentBearingSingle2| adds a |SpeedDependentBearing| configuration
% to the input parameter structure. The bearing stiffness and damping matrices
% are defined as [bearings x speed-points] lookup tables;
% |updateSpeedDependentMatrices| interpolates these tables at each frequency
% step during FRF analysis. The function:
% * Defines full 2x2 stiffness (Kxx, Kyy, Kxy, Kyx) and damping (Cxx, Cyy, Cxy, Cyx) matrices
% * Validates input data via |checkInputData|
% * Sorts bearings by axial position via |sortRowsWithShaftDis|
% * Sets |ComponentSwitch.hasSpeedDependentBearing = true|
%
%% Input Arguments
% * |InitialParameter| - Preconfigured system parameter structure created by an
%   |inputEssentialParameter*| function
%
%% Output Arguments
% * |Parameter| - Updated parameter structure with |SpeedDependentBearing| field appended
%
%% SpeedDependentBearing Parameters (SpeedDependentBearing structure)
% * |amount|                     - Number of speed-dependent bearings (scalar)
% * |inShaftNo|                  - Shaft index for each bearing [n×1 column vector]
% * |dofOfEachNodes|             - DOF per node; must be 0 (no bearing mass) [n×1 column vector]
% * |positionOnShaftDistance|    - Axial position from shaft left end [m] [n×1 column vector]
% * |speed|                      - Speed lookup vector [rad/s] [1×m row vector]
% * |stiffness|                  - Horizontal stiffness Kxx [N/m] [n×m matrix]
% * |stiffnessVertical|          - Vertical stiffness Kyy [N/m] [n×m matrix]
% * |stiffnessHV|                - Cross-coupled stiffness Kxy [N/m] [n×m matrix]
% * |stiffnessVH|                - Cross-coupled stiffness Kyx [N/m] [n×m matrix]
% * |damping|                    - Horizontal damping Cxx [Ns/m] [n×m matrix]
% * |dampingVertical|            - Vertical damping Cyy [Ns/m] [n×m matrix]
% * |dampingHV|                  - Cross-coupled damping Cxy [Ns/m] [n×m matrix]
% * |dampingVH|                  - Cross-coupled damping Cyx [Ns/m] [n×m matrix]
% * |mass|                       - Bearing mass [kg]; must be zero [n×1 column vector]
% * |isHertzian|                 - Hertzian contact flag; must be 0 for this type [n×1 column vector]
%
%% System Flags
% Automatically sets |ComponentSwitch.hasSpeedDependentBearing = true|.
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function Parameter = inputSpeedDependentBearingSingle2(InitialParameter)
% Input parameters for speed-dependent bearings
% (Linearized results considering non-linear factors like Hertzian contact)
SpeedDependentBearing.amount          = 2;
SpeedDependentBearing.inShaftNo       = [1; 1];
SpeedDependentBearing.dofOfEachNodes  = [0; 0]; % dof must be 0
SpeedDependentBearing.positionOnShaftDistance = [0; 517.2] * 10^-3; % from the left end of the shaft (m)

% Speed vector (rad/s)
SpeedDependentBearing.speed           = [100, 200, 300];

% Stiffness matrices (Rows: number of bearings, Columns: number of speed points)
SpeedDependentBearing.stiffness       = [1e6, 5e6, 1e7; 2e6, 7e6, 2e7]; % Horizontal stiffness Kxx (N/m)
SpeedDependentBearing.stiffnessVertical = [1e6, 5e6, 1e7; 2e6, 7e6, 2e7]; % Vertical stiffness Kyy (N/m)
SpeedDependentBearing.stiffnessHV     = [1e5, 5e5, 1e6; 2e5, 7e5, 2e6]; % Cross-coupled stiffness Kxy (N/m)
SpeedDependentBearing.stiffnessVH     = [1e5, 5e5, 1e6; 2e5, 7e5, 2e6]; % Cross-coupled stiffness Kyx (N/m)

% Damping matrices (Rows: number of bearings, Columns: number of speed points)
SpeedDependentBearing.damping         = [1e3, 2e3, 3e3; 4e3, 5e3, 6e3]; % Horizontal damping Cxx (N*s/m)
SpeedDependentBearing.dampingVertical = [1e3, 2e3, 3e3; 4e3, 5e3, 6e3]; % Vertical damping Cyy (N*s/m)
SpeedDependentBearing.dampingHV       = [1e2, 2e2, 3e2; 4e2, 5e2, 6e2]; % Cross-coupled damping Cxy (N*s/m)
SpeedDependentBearing.dampingVH       = [1e2, 2e2, 3e2; 4e2, 5e2, 6e2]; % Cross-coupled damping Cyx (N*s/m)

% Auxiliary parameters
SpeedDependentBearing.mass = zeros(SpeedDependentBearing.amount, 1); % kg
SpeedDependentBearing.isHertzian = zeros(SpeedDependentBearing.amount, 1); % boolean, must be 0

% Data validation and sorting
trans = rmfield(SpeedDependentBearing, 'speed');
checkInputData(trans)
trans = sortRowsWithShaftDis(trans);
trans.speed = SpeedDependentBearing.speed;

% Assign to main parameter struct
Parameter = InitialParameter;
Parameter.SpeedDependentBearing = trans;

% enable the flag
Parameter.ComponentSwitch.hasSpeedDependentBearing = true;
end
