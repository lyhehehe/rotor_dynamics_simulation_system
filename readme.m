%% SpoolDyn — Interactive Tutorial Script
%
% This script walks through the key features of SpoolDyn step by step.
% Run each section independently using MATLAB's "Run Section" feature
% (place the cursor inside a section and press Ctrl+Enter, or click
% "Run Section" in the Editor tab).
%
% Sections overview:
%   Section 1 - Single-shaft time-domain simulation (Example 1)
%   Section 2 - Twin-spool time-domain with Hertzian contact (Example 2)
%   Section 3 - Campbell diagram and critical speeds (Example 5)
%   Section 4 - Unbalance response frequency sweep (Example 6)
%   Section 5 - Speed-dependent bearings in FRF solver (Example 7)
%
% For full documentation see README.md.
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.

%% Section 1 — Single-shaft time-domain simulation
% Demonstrates the complete time-domain workflow for a simple two-disk
% single-shaft rotor (corresponding to main_example_1.m).
%
% Pipeline:
%   input -> establishModel -> generateDynamicEquation -> calculateResponse -> signalProcessing

clc; clear; close all;

% --- Step 1: Load input parameters ---
% Open inputEssentialParameterSingle2.m to view or edit shaft, disk,
% bearing, and operating-speed parameters.
InitialParameter = inputEssentialParameterSingle2();

% --- Step 2: Assemble FEM model ---
% establishModel() builds global matrices, creates mesh & model diagrams
% in the <meshDiagram> and <modelDiagram> folders.
Parameter = establishModel(InitialParameter);

% --- Step 3: Generate dynamic equation file ---
% Writes dynamicEquation.m to the project root. This file is tailored to
% the current model (DOF count, force terms, etc.) for efficient ODE solving.
generateDynamicEquation(Parameter);

% --- Step 4: Integrate in time ---
% Simulates 0 to 50 seconds at 2^14 samples/s using MATLAB's stiff solver.
% Returns displacement q [dofNum x nTimePoints], velocity dq, and time t.
[q, dq, t] = calculateResponse(Parameter, [0, 50], 2^14, 'calculateMethod', 'ode15s');

% --- Step 5: Post-process ---
% Specify which plots to generate (true = generate, false = skip).
% All figures are saved to the <signalProcess> folder.
SwitchFigure.displacement     = true;   % time history of X and Y displacements
SwitchFigure.axisTrajectory   = true;   % 2-D rotor orbit (X vs Y)
SwitchFigure.axisTrajectory3d = false;  % 3-D orbit
SwitchFigure.phase            = false;  % phase portrait
SwitchFigure.fftSteady        = true;   % steady-state FFT spectrum
SwitchFigure.fftTransient     = false;  % short-time Fourier transform (waterfall)
SwitchFigure.poincare         = false;  % Poincare section
SwitchFigure.poincare_phase   = false;  % Poincare + phase portrait
SwitchFigure.saveFig          = false;  % save .fig files (PNG always saved)
SwitchFigure.saveEps          = false;  % save .eps vector files

% Analyze only the steady-state window [25, 35] s to avoid run-up transients.
signalProcessing(q, dq, t, Parameter, [25, 35], 2^14, SwitchFigure)


%% Section 2 — Twin-spool time-domain with Hertzian contact
% Demonstrates a dual-shaft system with inter-shaft bearing and Hertzian
% rolling-element contact (corresponding to main_example_2.m).
%
% Key additions vs. Section 1:
%   - inputBearingHertzTwinSpool: defines Hertzian contact parameters
%   - inputIntermediateBearingTwinSpool: adds the inter-shaft bearing
%   - Manual mesh: user controls node density between key points
%   - ode15s is essential here because Hertzian contact creates stiff ODEs

clc; clear; close all;

% Load all three parameter files
InitialParameter = inputEssentialParameterTwinSpool();
InitialParameter = inputBearingHertzTwinSpool(InitialParameter);
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter);

% Manual mesh: each value is the number of sub-segments between two
% consecutive key nodes (disk/bearing positions) along the shaft.
manualGrid{1} = [1, 2, 1, 7, 1, 1, 3]; % shaft 1 — 7 key-node intervals
manualGrid{2} = [1, 3, 4, 3];           % shaft 2 — 4 key-node intervals
Parameter = establishModel(InitialParameter, 'gridFineness', manualGrid);

generateDynamicEquation(Parameter);

% Simulate 0-10 s. ode15s handles the stiff Hertzian contact forces.
[q, dq, t] = calculateResponse(Parameter, [0, 10], 2^14, 'calculateMethod', 'ode15s');

SwitchFigure.displacement     = true;
SwitchFigure.axisTrajectory   = false;
SwitchFigure.axisTrajectory3d = false;
SwitchFigure.phase            = false;
SwitchFigure.fftSteady        = false;
SwitchFigure.fftTransient     = false;
SwitchFigure.poincare         = false;
SwitchFigure.poincare_phase   = false;
SwitchFigure.saveFig          = false;
SwitchFigure.saveEps          = false;

signalProcessing(q, dq, t, Parameter, [0, 10], 2^14, SwitchFigure)


%% Section 3 — Campbell diagram and mode shapes
% Frequency-domain modal analysis: track natural frequencies across a
% speed sweep and visualize deformed rotor geometry at critical speeds.
% Corresponds to main_example_5_campbell_and_mode.m.
%
% This path does NOT require generateDynamicEquation or calculateResponse.

clc; clear; close all;

% Build the twin-spool model (same setup as Section 2)
InitialParameter = inputEssentialParameterTwinSpool();
InitialParameter = inputBearingHertzTwinSpool(InitialParameter);
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter);

manualGrid{1} = [1, 2, 1, 7, 1, 1, 3];
manualGrid{2} = [1, 3, 4, 3];
Parameter = establishModel(InitialParameter, 'gridFineness', manualGrid, ...
            'isPlotMesh', false, 'isPlotModel', false); % suppress auto-plots

% --- Campbell Diagram ---
% exciteRad: the shaft-1 speed sweep (rad/s). Natural frequencies are
% computed at each point, and crossings with the 1x excitation line
% identify critical speeds.
max_rpm   = 10000;
exciteRad = linspace(1, max_rpm/60*2*pi, 500);

[eigMatrix, criticalSpeeds] = calculateCampbell(Parameter, exciteRad, ...
    'isPlot',          true,    ... % draw Campbell diagram with critical-speed markers
    'isFilter',        true,    ... % track mode lines to avoid artefact crossings
    'filterMethod',    'slope', ... % 'slope' (gradient) or 'MAC' (correlation)
    'isUseGyroMatrix', true);       % include gyroscopic effect (recommended)

% eigMatrix   — [nModes x 500] natural frequencies in rad/s
% criticalSpeeds — sorted 1× crossing speeds [rad/s]
fprintf('Found %d critical speeds (rad/s): ', length(criticalSpeeds));
fprintf('%.1f  ', criticalSpeeds); fprintf('\n');

% --- Mode Shapes at Critical Speeds ---
% Renders the deformed shaft geometry at each critical speed.
% scaleFactor = 'auto' scales max deflection to 1.5x the shaft outer radius.
[ModeShapes, ZCoords] = calculateModeShape(Parameter, criticalSpeeds, ...
    'isPlot',          true,  ...
    'direction',       'X',   ... % 'X' (horizontal) or 'Y' (vertical)
    'isUseGyroMatrix', true);

% ModeShapes{iShaft, iSpeed} contains the displacement at each shaft node.
% ZCoords{iShaft} contains the global Z-axis coordinates of those nodes.


%% Section 4 — Unbalance response frequency sweep
% Solves the synchronous FRF at each operating condition and returns the
% complex steady-state displacement vector. Corresponds to
% main_example_6_solve_in_frequency_domain.m.

clc; clear; close all;

% Build model (suppress diagrams for speed)
InitialParameter = inputEssentialParameterTwinSpool();
InitialParameter = inputBearingHertzTwinSpool(InitialParameter);
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter);

manualGrid{1} = [1, 2, 1, 7, 1, 1, 3];
manualGrid{2} = [1, 3, 4, 3];
Parameter = establishModel(InitialParameter, 'gridFineness', manualGrid, ...
            'isPlotMesh', false, 'isPlotModel', false);

% --- Usage 1: Single condition from Parameter.Status ---
% Returns response [dofNum x 1 x shaftNum].
% response(:, 1, 1) — vibration excited by shaft-1 unbalance at shaft-1's speed.
% response(:, 1, 2) — vibration excited by shaft-2 unbalance at shaft-2's speed.
response1 = calculateUnbalanceResponse(Parameter);

% --- Usage 2: Multi-condition speed sweep ---
% speedMatrix rows = shafts, columns = operating conditions (rad/s).
% Returns [dofNum x numConditions x shaftNum].
speedMatrix = [100, 200, 300;   % shaft 1 speeds at conditions 1, 2, 3
               200, 400, 600];  % shaft 2 speeds
response2 = calculateUnbalanceResponse(Parameter, speedMatrix);

% Example: extract amplitude and phase at DOF 1, shaft-1 excitation
dof = 1;
amplitude_mm = abs(response2(dof, :, 1)) * 1e3; % convert m → mm
phase_deg    = angle(response2(dof, :, 1)) * 180/pi;
fprintf('Node 1 (shaft-1 excitation): amplitude at 3 conditions = %.3e, %.3e, %.3e mm\n', ...
        amplitude_mm(1), amplitude_mm(2), amplitude_mm(3));

% --- Usage 3: Sweep with orbit trajectory visualization ---
% Superimposes contributions from both shafts into a physical orbit.
speedMatrix_wide = [linspace(10, 500, 100);   % shaft 1
                    linspace(20, 1000, 100)];  % shaft 2 (ratio 2:1)
response3 = calculateUnbalanceResponse(Parameter, speedMatrix_wide, ...
    'isPlot',           true, ...
    'plotNodeID',       [1, 4], ... % node IDs to visualize
    'plotConditionIdx', [50, 100]); % which columns of speedMatrix to plot


%% Section 5 — Speed-dependent bearings in FRF solver
% Demonstrates bearings whose stiffness and damping coefficients are
% tabulated as functions of speed (e.g., from linearized Hertzian contact
% analysis). Corresponds to main_example_7_speed_independent_bearing.m.
%
% How it works:
%   1. inputSpeedDependentBearingSingle2 stores lookup tables in Parameter.SpeedDependentBearing
%   2. The standard Bearing struct is zeroed out so it does not contribute
%   3. calculateUnbalanceResponse calls updateSpeedDependentMatrices at each
%      frequency step to interpolate K and C from the lookup tables

clc; clear; close all;

% Load essential parameters for a single-shaft system
InitialParameter = inputEssentialParameterSingle2();

% Add the speed-dependent bearing configuration.
% Open inputSpeedDependentBearingSingle2.m to edit the speed lookup
% vectors and the [n_bearings x n_speeds] K/C matrices.
InitialParameter = inputSpeedDependentBearingSingle2(InitialParameter);

% Remove the standard (constant) bearings — they are replaced by the
% speed-dependent ones defined above.
InitialParameter.Bearing.amount                   = 0;
InitialParameter.Bearing.inShaftNo                = [];
InitialParameter.Bearing.dofOfEachNodes           = [];
InitialParameter.Bearing.positionOnShaftDistance  = [];
InitialParameter.Bearing.stiffness                = [];
InitialParameter.Bearing.stiffnessVertical        = [];
InitialParameter.Bearing.damping                  = [];
InitialParameter.Bearing.dampingVertical          = [];
InitialParameter.Bearing.mass                     = [];
InitialParameter.Bearing.isHertzian               = [];

% Assemble the model
Parameter = establishModel(InitialParameter);

% Compute unbalance response across a speed sweep.
% At each of the 100 speed conditions, calculateUnbalanceResponse calls
% updateSpeedDependentMatrices to interpolate the bearing K and C tables
% at the current operating speed before solving H(omega)*X = F.
speedMatrix = linspace(0, 500, 100); % single shaft: 1 x 100 row vector

response = calculateUnbalanceResponse(Parameter, speedMatrix, ...
    'isPlot',           true,      ...
    'plotNodeID',       [1, 4],    ... % plot orbits at nodes 1 and 4
    'plotConditionIdx', [50, 100]);    % at the 50th and 100th speed condition

% Plot a simple amplitude vs. speed Bode curve for node 1 (DOF 1)
figure('Name', 'Amplitude vs. Speed — Node 1');
nodeStartDof = Parameter.Mesh.dofInterval(1, 1);  % first DOF of node 1
speed_rpm = speedMatrix * 30/pi;                  % rad/s → RPM
amp_um    = abs(response(nodeStartDof, :, 1)) * 1e6; % m → µm

plot(speed_rpm, amp_um, 'LineWidth', 1.5);
xlabel('Speed (RPM)');
ylabel('Amplitude (\mum)');
title('Unbalance Response at Node 1 — Speed-Dependent Bearing');
grid on;
