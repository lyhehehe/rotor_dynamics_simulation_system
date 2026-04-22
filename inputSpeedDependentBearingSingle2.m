%% inputIntermediateBearing - Configure intermediate bearings with Hertzian contact
%
% This function configures parameters for intermediate bearings connecting 
% dual-shaft systems, supporting both linear and Hertzian contact models 
% for rotor dynamics analysis.
%
%% Syntax
%  OutputParameter = inputIntermediateBearing(InputParameter)
%
%% Description
% |inputIntermediateBearing| adds intermediate bearing configuration to 
% existing rotor system parameters. It supports:
% * Linear spring-damper connections between shafts
% * Mass-spring chain connections
% * Hertzian contact force modeling
%
% * Inputs:
%   * |InputParameter| - Preconfigured system parameters structure 
%
% * Outputs:
%   * |OutputParameter| - Updated parameter structure with intermediate bearings
%
%% Intermediate Bearing Parameters (IntermediateBearing structure)
% * amount              - Number of intermediate bearings (scalar)
% * betweenShaftNo      - Connected shaft indices [n×2 matrix]
% * dofOfEachNodes      - Degrees of freedom per node (column vector)
% * positionOnShaftDistance - Mounting positions from shaft ends [n×2 matrix, m]
% * isHertzian          - Hertzian contact activation flags (logical column)
% * isHertzianTop       - Hertzian force position flags (logical column)
% * stiffness           - Horizontal stiffness [N/m] (matrix)
% * stiffnessVertical   - Vertical stiffness [N/m] (matrix)
% * damping             - Horizontal damping [Ns/m] (matrix)
% * dampingVertical     - Vertical damping [Ns/m] (matrix)
% * mass                - Intermediate masses [kg] (column vector)
% * rollerNum           - Number of rolling elements (column vector)
% * radiusInnerRace     - Inner race radii [m] (column vector)
% * radiusOuterRace     - Outer race radii [m] (column vector)
% * innerShaftNo        - Shaft containing inner race (column vector)
% * clearance           - Bearing clearances [m] (column vector)
% * contactStiffness    - Hertzian stiffness [N/m^1.5] (column vector)
% * coefficient         - Contact force exponent (column vector)
%
%% Model Configuration Rules
% * Connection Types:
%   * Basic Connection: Linear spring-damper between shafts
%   * Mass-spring Chain: Multiple masses with sequential stiffness/damping
%   * Hertzian Contact: Nonlinear force at specified position 
%     (isHertzianTop=true for top connection, false for bottom)
% * Automatic Sorting:
%   * Shaft indices automatically sorted in ascending order
%   * Associated parameters (stiffness, mass, DOF) reordered consistently
%
%% System Flags
% Automatically enables:
% * |hasIntermediateBearing| in ComponentSwitch
% * |hasHertzianForce| if any bearing has |isHertzian=true|
%
%% Example
%   % Initialize system parameters
%   sysParams = inputEssentialParameterBO();
%   % Configure intermediate bearings
%   sysParams = inputIntermediateBearing(sysParams);
%
%% See Also
%  checkInputData, sortRowsWithShaftDis, inputEssentialParameterBO
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


%%
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

