%% inputEssentialParameterTwinSpool2 - Input essential parameters for a dual-rotor system
%
% This function initializes and returns the essential parameters for 
% modeling a rotor system, including shaft configurations, disk 
% properties, system status, and component switches. The parameters are 
% organized in a structured format for dynamic analysis.
%
%% Syntax
%  InitialParameter = inputEssentialParameterTwinSpool2()
%
%% Description
% |inputEssentialParameterTwinSpool2()| initializes parameters for shafts, disks, 
% system running status, and component configuration switches. These 
% parameters are essential for building the mathematical model of a 
% rotor system in subsequent dynamic analysis.
%
%% Output Parameters
% * InitialParameter - Structure containing all system parameters with fields:
%   * Status       - System running status parameters (acceleration, speed profile)
%   * Shaft        - Geometric and material properties of shafts
%   * Disk         - Geometric and inertial properties of disks
%   * ComponentSwitch - Boolean flags for system component activation
%
%% Shaft Parameters (Shaft structure)
% * amount           - Number of shafts (scalar)
% * totalLength      - Column vector of shaft lengths [m]
% * dofOfEachNodes   - Degrees of freedom per node (column vector)
% * segmentLength    - Segment lengths per shaft [cell array of vectors]
% * outerRadius      - Physical outer radii [m] [cell array of vectors]
% * innerRadius      - Physical inner radii [m] [cell array of vectors]
% * outerRadiusStiff - Stiffness outer radii [m] [cell array of vectors]
% * innerRadiusStiff - Stiffness inner radii [m] [cell array of vectors]
% * density          - Material densities [kg/m³] [cell array of vectors]
% * elasticModulus   - Elastic moduli [Pa] [cell array of vectors]
% * poissonRatio     - Poisson's ratios [cell array of vectors]
% * eccentricity     - Distributed mass eccentricity [m] [cell array of vectors]
% * eccentricityPhase - Phase of distributed eccentricity [rad] [cell array of vectors]
% * rayleighDamping  - Rayleigh damping coefficients [alpha, beta]
%
%% Disk Parameters (Disk structure)
% * amount               - Number of disks (scalar)
% * inShaftNo            - Shaft index for each disk (column vector)
% * dofOfEachNodes       - DOF per node (column vector)
% * innerRadius          - Disk inner radii [m] (column vector)
% * outerRadius          - Disk outer radii [m] (column vector)
% * thickness            - Disk thicknesses [m] (column vector)
% * positionOnShaftDistance - Mounting positions from shaft ends [m] (column)
% * density              - Material densities [kg/m³] (column vector)
% * eccentricity         - Mass eccentricities [m] (column vector)
% * eccentricityPhase    - Phase of mass eccentricities [rad] (column vector)
%
%% Status Parameters (Status structure)
% * ratio            - Speed ratio between shafts (vector)
% * vmax             - Maximum rotational speed of shaft 1 [rad/s]
% * acceleration     - Rotational acceleration of shaft 1 [rad/s²]
% * duration         - Duration at maximum speed [s]
% * isDeceleration   - Flag for deceleration phase (logical)
% * vmin             - Minimum speed after deceleration [rad/s]
% * isUseCustomize   - Flag for custom speed profile (logical)
% * customize        - Handle to custom speed profile function
%
%% Component Switches (ComponentSwitch structure)
% * hasGravity              - Enable gravitational effects (logical)
% * hasIntermediateBearing  - Enable intermediate bearings (logical)
% * hasLoosingBearing       - Enable bearing clearance (logical)
% * hasRubImpact            - Enable rotor-stator rub (logical)
% * hasCouplingMisalignment - Enable coupling misalignment (logical)
% * hasHertzianForce        - Enable Hertzian contact forces (logical)
% * hasCustom               - Enable custom component (logical)
%
%% Example
%   InitialParams = inputEssentialParameterTwinSpool();
%   % Access shaft parameters:
%   shaftLengths = InitialParams.Shaft.totalLength;
%
%% See Also
%  checkInputData, calculateStatus
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


%%
function InitialParameter = inputEssentialParameterTwinSpool2()

% typing the parameter about shaft
Shaft.amount            = 2;
Shaft.totalLength       = [962; 382]*10^-3; % all vectors in column (m)
Shaft.dofOfEachNodes    = 4 * ones(Shaft.amount,1);
Shaft.segmentLength     = {[230, 10, 10, 10, 10, 10, 682]*10^-3;
                           [270, 10, 10, 10, 10, 10, 62]*10^-3}; % n*1 cell saving segment Length of each shaft 
Shaft.outerRadius       = {[10,   20, 40, 80, 40, 20, 10]*10^-3;
                           [32.5, 40, 60, 80, 60, 40, 32.5]*10^-3}; % m
Shaft.innerRadius       = {[0,  0,  10, 20, 10, 0,  0]*10^-3;
                           [20, 20, 20, 20, 20, 20, 20]*10^-3}; % m
Shaft.outerRadiusStiff  = {[10,   15, 20, 40, 20, 15, 10]*10^-3;
                           [32.5, 35, 50, 60, 50, 35, 32.5]*10^-3}; % m
Shaft.innerRadiusStiff  = {[0,  0,  0, 20, 0, 0,  0]*10^-3;
                           [20, 20, 30, 40, 30, 20, 20]*10^-3}; % m
Shaft.density           = {7850 * ones(1,7); 7850 * ones(1,7)}; % kg/m^3
Shaft.elasticModulus    = {210e9 * ones(1,7); 210e9 * ones(1,7)}; % Pa
Shaft.poissonRatio      = {0.296 * ones(1,7); 0.296 * ones(1,7)};
Shaft.eccentricity      = {[0, 0, 0, 0.0979e-3, 0, 0, 0]; [0, 0, 0, 0.0979e-3, 0, 0, 0]}; % m
Shaft.eccentricityPhase = {zeros(1,7); zeros(1,7)}; % rad
checkInputData(Shaft)
Shaft.rayleighDamping   = [0, 3e-4]; % [alpha, beta] CShaft = alpha*(MShaft+MDisk) + beta*KShaft

%%

% typing the parameter about running status
Status.ratio            = [1.3]; % [v-shaft2/v-shaft1; v-shaft3/v-shaft1]
Status.vmax             = 200; % rad/s, the maximum rotational speed for shaft 1
Status.acceleration     = 20; % rad/s^2, acceleration of shaft 1
Status.duration         = 0; % s, the duration of shaft 1 in vmax
Status.isDeceleration   = true; % boolean, add a deceleration in status
Status.vmin             = 0; % s, the minimum speed afterdeceleration

% (otherwise) you can define your own simulation status function
% define your own function in calculateStatus() where the single time point
% is input, the output must be [acceleration, speed, angular] corresponding
% to input time "tn", in each output the dimension is m*X, X is the number
% of the shaft, m is number of elements in tn
Status.isUseCustomize   = false;
Status.customize        = @(tn) calculateStatus(tn);

% check input
if  length(Status.ratio) >= Shaft.amount
    error('too much input parameter in Status.ratio')
end

%%

% typing the parameter about disk
Disk.amount             = 2;
Disk.inShaftNo          = [1, 2]'; % disks in the i-th shaft
Disk.dofOfEachNodes     = 4 * ones(Disk.amount,1);
Disk.innerRadius        = [10, 32.5]' *10^-3; % m
Disk.outerRadius        = [125, 125]' *10^-3; % m
Disk.thickness          = 0.015*ones(Disk.amount,1); % m
Disk.positionOnShaftDistance = [718.5, 150.5]' * 10^-3; %from left end (m)
Disk.density            = 7850*ones(Disk.amount,1); % kg/m^3
Disk.eccentricity       = 0.0979e-3*ones(Disk.amount,1); % m
Disk.eccentricityPhase  = zeros(Disk.amount, 1); % rad/s

% check input
checkInputData(Disk)

for iDisk = 1:1:Disk.amount
   if Shaft.dofOfEachNodes(Disk.inShaftNo(iDisk)) ~= Disk.dofOfEachNodes(iDisk)
      error(['the dof of each disk should equal to the dof of the shaft'...
            ,' this disk locating']); 
   end
end

%%

% typing the parameter about linear bearing
% If you choose to input the bearing parameter here, you should not use the
% inputBearingHertz()
% model: shaft--k1c1--mass--k2c2--basement

Bearing.amount          = 0;
Bearing.inShaftNo       = [];
Bearing.dofOfEachNodes  = []; % if mass=0, dof must be 0 
Bearing.positionOnShaftDistance = []; % m
Bearing.stiffness       = []; % N*m
Bearing.stiffnessVertical = []; % N*m
Bearing.damping         = []; % N*s/m
Bearing.dampingVertical = []; % N*s/m
Bearing.mass            = []; % kg
Bearing.isHertzian      = [];

checkInputData(Bearing)

%%

% ComponentSwitch will be changed by corresponding input..() function
ComponentSwitch.hasGravity = true; % boolean, true-> take gravity into account; false->no gravity in the dynamic equation
ComponentSwitch.hasIntermediateBearing = false;
ComponentSwitch.hasLoosingBearing = false;
ComponentSwitch.hasRubImpact = false;
ComponentSwitch.hasCouplingMisalignment = false;
ComponentSwitch.hasHertzianForce = false;
ComponentSwitch.hasCustom = false;


%%

% Output initialParameter without optional parameter
InitialParameter.Status          = Status;
InitialParameter.Shaft           = Shaft;
InitialParameter.Disk            = Disk;
InitialParameter.Bearing         = Bearing;
InitialParameter.ComponentSwitch = ComponentSwitch;


end