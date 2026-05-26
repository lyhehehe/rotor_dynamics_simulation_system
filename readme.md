# SpoolDyn — Multispool Rotordynamics Simulator

**SpoolDyn** is an open-source MATLAB tool for simulating, analyzing, and post-processing the vibration of rotating machinery — from a simple Jeffcott rotor to complex twin-spool aero-engine configurations. It is built on the Finite Element Method (Timoshenko beams), lumped-parameter modeling, and automatic code generation for fast nonlinear time integration.

> **Requires MATLAB R2025a or later.**

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Workflow Overview](#4-workflow-overview)
5. [Examples](#5-examples)
6. [User Guide](#6-user-guide)
   - [6.1 Input Module](#61-input-module)
   - [6.2 Modeling Module](#62-modeling-module)
   - [6.3 Time-Domain Analysis](#63-time-domain-analysis)
   - [6.4 Frequency-Domain Analysis](#64-frequency-domain-analysis)
   - [6.5 Post-Processing](#65-post-processing)
7. [References](#7-references)

---

## 1. Overview

SpoolDyn uses four fundamental element types — **shaft**, **disk**, **bearing**, and **inter-shaft bearing** — assembled into a global FEM model. Users interact through a small set of function calls; SpoolDyn handles matrix assembly, code generation, numerical integration, and visualization automatically.

**Key capabilities:**

| Feature | Description |
|---|---|
| Multi-shaft systems | Unlimited co-axial or coupled shafts |
| Timoshenko beam theory | Shear deformation and rotary inertia included |
| Gyroscopic effects | Speed-dependent gyroscopic matrix, fully coupled |
| Nonlinear Hertzian contact | Rolling-element bearing forces, ball and roller bearings |
| Inter-shaft bearings | Connect any two shafts at arbitrary axial positions |
| Speed-dependent bearings | Linearized stiffness/damping lookup tables, interpolated per FRF step |
| Campbell diagram | Tracks natural frequencies vs. speed, marks critical speeds |
| Mode shapes | Deformed geometry visualization at critical speeds |
| Unbalance response (FRF) | Steady-state synchronous response across a speed sweep |
| Time-domain integration | RK4, `ode45`, `ode15s`, `ode23s`; run-up / run-down profiles |
| Auto code generation | Generates `dynamicEquation.m` tailored to the model for fast ODE solving |

---

## 2. Installation

### Requirements

- **MATLAB R2025a or later** (no additional toolboxes required)

### Method 1 — Git (Recommended)

```bash
git clone https://github.com/AlkaidWood/spooldyn.git
cd spooldyn
```

Then add the folder to your MATLAB path:

```matlab
addpath(genpath('path/to/spooldyn'))
```

### Method 2 — ZIP Download

1. Go to [https://github.com/AlkaidWood/spooldyn](https://github.com/AlkaidWood/spooldyn)
2. Click **Code → Download ZIP** and extract to any folder.
3. Add the folder to your MATLAB path.

---

## 3. Quick Start

Run one of the ready-made example scripts from the MATLAB command window:

```matlab
cd path/to/spooldyn
run('main_example_1.m')   % single-shaft, time-domain run-up
```

Output diagrams are saved automatically in `modelDiagram/`, `meshDiagram/`, and `signalProcess/`.

For a step-by-step walkthrough of every major feature, open the companion tutorial:

```matlab
open('readme.m')          % run each section with Ctrl+Enter
```

---

## 4. Workflow Overview

SpoolDyn supports two main analysis paths:

### Path A — Time-Domain Analysis

```
inputEssentialParameter*()   →  Define geometry, material, operating speed
establishModel()             →  FEM assembly, mesh, global matrices
generateDynamicEquation()    →  Code-generate ODE function file
calculateResponse()          →  Numerical time integration (run-up, steady state, run-down)
signalProcessing()           →  FFT, orbits, Poincaré, STFT, save figures
```

Best for: transient run-ups, nonlinear Hertzian contact, rub-impact, custom external forces.

### Path B — Frequency-Domain Analysis

```
inputEssentialParameter*()   →  Define geometry, material, operating speed
establishModel()             →  FEM assembly, mesh, global matrices
calculateCampbell()          →  Natural frequencies vs. speed → critical speeds
calculateModeShape()         →  Deformed geometry at critical speeds
calculateUnbalanceResponse() →  Synchronous FRF sweep (Bode / waterfall)
```

Best for: design-phase screening, critical speed prediction, unbalance sensitivity.

> **You do not need `generateDynamicEquation` or `calculateResponse` for Path B.**

---

## 5. Examples

Seven ready-to-run scripts are provided in the project root:

| Script | System | Analysis | New features |
|---|---|---|---|
| `main_example_1.m` | Two-disk, single shaft | Time-domain | Basic workflow |
| `main_example_2.m` | Twin-spool, Hertzian contact | Time-domain | Multi-shaft, inter-shaft bearing, Hertzian |
| `main_example_3.m` | Variable cross-section shaft | Time-domain | Cross-coupled bearings |
| `main_example_4.m` | Twin-spool, inter-shaft bearing | Time-domain | Intermediate bearing |
| `main_example_5_campbell_and_mode.m` | Twin-spool | Frequency-domain | **Campbell diagram, mode shapes** |
| `main_example_6_solve_in_frequency_domain.m` | Twin-spool | Frequency-domain | **Unbalance response sweep, orbit plots** |
| `main_example_7_speed_independent_bearing.m` | Single shaft | Frequency-domain | **Speed-dependent bearing** |


### Example 1 — Single-Shaft Run-Up

```matlab
InitialParameter = inputEssentialParameterSingle2();

Parameter = establishModel(InitialParameter);
generateDynamicEquation(Parameter);
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'calculateMethod', 'ode15s');

SwitchFigure.displacement     = true;
SwitchFigure.axisTrajectory   = false;
SwitchFigure.axisTrajectory3d = false;
SwitchFigure.phase            = false;
SwitchFigure.fftSteady        = false;
SwitchFigure.fftTransient     = false;
SwitchFigure.poincare         = false;
SwitchFigure.poincare_phase   = false;
SwitchFigure.saveFig          = true;
SwitchFigure.saveEps          = false;
signalProcessing(q, dq, t, Parameter, [0,50], 2^14, SwitchFigure)
```

![example1-mesh](readme/example1/MeshResultOfShaft1.png)

![example1-model](readme/example1/theWholeModel.png)

### Example 2 — Twin-Spool with Hertzian Contact

```matlab
InitialParameter = inputEssentialParameterTwinSpool();
InitialParameter = inputBearingHertzTwinSpool(InitialParameter);
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter);

manualGrid{1} = [1,2,1,7,1,1,3];  % mesh for shaft 1
manualGrid{2} = [1,3,4,3];         % mesh for shaft 2
Parameter = establishModel(InitialParameter, 'gridFineness', manualGrid);
generateDynamicEquation(Parameter);
[q, dq, t] = calculateResponse(Parameter, [0,10], 2^14, 'calculateMethod', 'ode15s');
signalProcessing(q, dq, t, Parameter, [0,10], 2^14, SwitchFigure)
```

![example2-mesh1](readme/example2/MeshResultOfShaft1.png)
![example2-mesh2](readme/example2/MeshResultOfShaft2.png)

Model diagram of shaft 1
![example2-model1](readme/example2/diagramOfShaft1.png)
Model diagram of shaft 2
![example2-model2](readme/example2/diagramOfShaft2.png)
Full rotor model
![example2-model](readme/example2/theWholeModel.png)

---

## 6. User Guide

### 6.1 Input Module

All physical parameters are stored in dedicated MATLAB function files. SpoolDyn provides these input function types:

| Function | Required? | Purpose |
|---|---|---|
| `inputEssentialParameter*()` | **Yes** | Shafts, disks, basic bearings, operating speed |
| `inputBearingHertz*()` | Optional | Hertzian rolling-element bearing parameters |
| `inputIntermediateBearing*()` | Optional | Inter-shaft (inter-spool) bearing parameters |
| `inputSpeedDependentBearing*()` | Optional | Speed-tabulated stiffness/damping for FRF analysis |
| `inputCustomFunction()` | Optional | Custom external force function |

Call them in sequence before `establishModel`:

```matlab
InitialParameter = inputEssentialParameter*();
InitialParameter = inputBearingHertz*(InitialParameter);        % optional
InitialParameter = inputIntermediateBearing*(InitialParameter); % optional
Parameter = establishModel(InitialParameter);
```

> **Note:** Define bearing stiffness/damping either in `inputEssentialParameter*` or in `inputBearingHertz*`, but **not both**.

---

#### Shafts

Single-shaft configuration:

```matlab
Shaft.amount         = 1;
Shaft.totalLength    = 517.2e-3;       % m
Shaft.dofOfEachNodes = 4;              % always 4 (x, θx, y, θy per node)
Shaft.outerRadius    = 5e-3;           % m
Shaft.innerRadius    = 0;              % m (0 = solid shaft)
Shaft.density        = 7850;           % kg/m³
Shaft.elasticModulus = 207e9;          % Pa
Shaft.poissonRatio   = 0.3;
checkInputData(Shaft)
Shaft.rayleighDamping = [0, 1.9e-4];  % [α, β]: C = α·M + β·K (rotating parts only)
```

Twin-shaft configuration — all parameters as column vectors, one row per shaft:

```matlab
Shaft.amount         = 2;
Shaft.totalLength    = [962; 382]*1e-3;
Shaft.outerRadius    = [10; 32.5]*1e-3;
Shaft.innerRadius    = [0; 20]*1e-3;
Shaft.density        = [7850; 7850];
Shaft.elasticModulus = [210e9; 210e9];
Shaft.poissonRatio   = [0.296; 0.296];
Shaft.dofOfEachNodes = [4; 4];
checkInputData(Shaft)
Shaft.rayleighDamping = [0, 3e-4];
```

> Shear modulus is auto-computed from `poissonRatio` and `elasticModulus`. SpoolDyn supports any number of shafts.

---

#### Disks

```matlab
Disk.amount                   = 4;
Disk.inShaftNo                = [1; 1; 2; 2];              % which shaft each disk is on
Disk.dofOfEachNodes           = 4*ones(4,1);
Disk.innerRadius              = [10; 10; 32.5; 32.5]*1e-3; % m
Disk.outerRadius              = [125; 125; 125; 125]*1e-3; % m
Disk.thickness                = 0.015*ones(4,1);           % m
Disk.positionOnShaftDistance  = [270.5; 718.5; 150.5; 292.5]*1e-3; % from shaft left end (m)
Disk.density                  = 7850*ones(4,1);            % kg/m³
Disk.eccentricity             = 0.0979e-3*ones(4,1);       % mass eccentricity (m)
Disk.eccentricityPhase        = zeros(4,1);                % phase (rad)
```

> `positionOnShaftDistance` is measured from the **left end** of the shaft each disk is mounted on.

---

#### Bearings

**Case 1 — Massless linear bearings** (defined inside `inputEssentialParameter*`):

```matlab
Bearing.amount                  = 2;
Bearing.inShaftNo               = [1; 1];
Bearing.dofOfEachNodes          = [0; 0];          % 0 = no bearing mass node
Bearing.positionOnShaftDistance = [0; 517.2]*1e-3; % m
Bearing.stiffness               = [1e7; 2e7];      % horizontal Kxx (N/m)
Bearing.stiffnessVertical       = [1e7; 2e7];      % vertical Kyy (N/m)
Bearing.damping                 = [1e4; 2e4];      % horizontal Cxx (Ns/m)
Bearing.dampingVertical         = [1e4; 2e4];      % vertical Cyy (Ns/m)
Bearing.mass                    = [0; 0];          % kg
Bearing.isHertzian              = [0; 0];
```

<img src="readme/bearing_model/massless-bearing.png" width="100%" height="" />

**Case 2 — Bearings with mass and Hertzian contact** (defined inside `inputBearingHertz*`):

```matlab
Bearing.amount                  = 2;
Bearing.inShaftNo               = [1; 1];
Bearing.dofOfEachNodes          = [2; 2];                    % 2 DOFs per mass node
Bearing.positionOnShaftDistance = [0; 517.2]*1e-3;
Bearing.stiffness               = [1e7, 1e8; 2e7, 2e8];     % [n × (mass_nodes+1)]
Bearing.stiffnessVertical       = [1e7, 1e8; 2e7, 2e8];
Bearing.damping                 = [1e4, 1e5; 2e4, 2e5];
Bearing.dampingVertical         = [1e4, 1e5; 2e4, 2e5];
Bearing.mass                    = [1; 1];                    % kg
Bearing.isHertzian              = [true; true];
Bearing.rollerNum               = [16; 16];
Bearing.radiusInnerRace         = [1e-2; 1e-2];              % m
Bearing.radiusOuterRace         = [1.5e-2; 1.5e-2];          % m
Bearing.clearance               = [1e-7; 1e-7];              % m
Bearing.contactStiffness        = [1e8; 1e8];                % N/m^(3/2) ball; N/m^(10/9) roller
Bearing.coefficient             = [3/2; 3/2];                % 3/2 = ball; 10/9 = roller
```

<img src="readme/bearing_model/mass-bearing.png" width="100%" height="" />

**Case 3 — Mixed bearings** (zero-pad to the maximum mass-block count):

```matlab
Bearing.amount        = 4;
Bearing.dofOfEachNodes = [0, 0, 0;   % bearing 1: no mass
                          2, 0, 0;   % bearing 2: 1 mass block
                          2, 2, 0;   % bearing 3: 2 mass blocks
                          2, 2, 2];  % bearing 4: 3 mass blocks
% Stiffness, damping, and mass matrices: zero-pad unused columns
```

<img src="readme/bearing_model/mixed-bearing.png" width="100%" height="" />

> **SpoolDyn automatically re-orders bearings and disks by axial position.**

---

#### Inter-shaft Bearings

Inter-shaft bearings connect two shafts at specified axial positions:

```matlab
IntermediateBearing.amount                   = 1;
IntermediateBearing.betweenShaftNo           = [1, 2];    % connects shaft 1 ↔ shaft 2
IntermediateBearing.dofOfEachNodes           = 2;         % 0 = massless, 2 = with mass
IntermediateBearing.positionOnShaftDistance  = [0.4, 0.2]; % [pos_on_shaft1, pos_on_shaft2] (m)
IntermediateBearing.isHertzian               = true;
IntermediateBearing.isHertzianTop            = true;
IntermediateBearing.stiffness                = [0, 1e7];  % N/m
IntermediateBearing.stiffnessVertical        = [0, 1e7];
IntermediateBearing.damping                  = [0, 0];
IntermediateBearing.dampingVertical          = [0, 0];
IntermediateBearing.mass                     = 0.3;       % kg
IntermediateBearing.innerShaftNo             = 1;         % which shaft is the inner rotor
IntermediateBearing.rollerNum                = 15;
IntermediateBearing.radiusInnerRace          = 15e-3;     % m
IntermediateBearing.radiusOuterRace          = 25e-3;     % m
IntermediateBearing.clearance                = 5e-6;      % m
IntermediateBearing.contactStiffness         = 13.34e9;   % N/m^(3/2)
IntermediateBearing.coefficient              = 3/2;
```

`betweenShaftNo` rows encode the connection topology:
`[1,2; 2,3; 1,3]` means: bearing 1 connects shafts 1–2, bearing 2 connects shafts 2–3, bearing 3 connects shafts 1–3.

<img src="readme/inter_shaft_bearing_model/inter-shaft-bearing.png" width="100%" height="" />

---

#### Operation Status

**Method 1 — Built-in speed profiles:**

```matlab
Status.ratio          = [1.3];   % [ω₂/ω₁; ω₃/ω₁; ...] — empty [] for single shaft
Status.vmax           = 200;     % rad/s — max speed of shaft 1
Status.acceleration   = 20;      % rad/s² — constant acceleration
Status.duration       = 0;       % s — dwell time at vmax (0 = no dwell)
Status.isDeceleration = true;    % ramp back down after vmax
Status.vmin           = 0;       % rad/s — final speed after deceleration
Status.isUseCustomize = false;
```

<img src="readme/rotation_status/rotation_status.png" width="100%" height="auto" />

**Method 2 — Custom speed function:**

```matlab
Status.isUseCustomize = true;
Status.customize = @(tn) calculateStatus(tn);  % returns [phase, speed, accel] per shaft
```

---

#### Custom External Forces

```matlab
Custom.amount                  = 2;
Custom.inShaftNo               = [1; 1];
Custom.positionOnShaftDistance = [20; 497.2]*1e-3; % m
Custom.force = @(qn, dqn, tn, omega, domega, ddomega, Parameter) ...
               customForce(qn, dqn, tn, omega, domega, ddomega, Parameter);
```

Edit `customForce.m` to return a `dofNum×1` force vector at each time step.

---

### 6.2 Modeling Module

```matlab
Parameter = establishModel(InitialParameter);
```

Outputs: `Parameter.Matrix` (sparse global matrices), `Parameter.Mesh.Node` (node data), `Parameter.Mesh.dofInterval` (DOF map).

#### Mesh Control

```matlab
% Automatic mesh
Parameter = establishModel(InitialParameter, 'gridFineness', 'low');    % key nodes only (default)
Parameter = establishModel(InitialParameter, 'gridFineness', 'middle'); % adds intermediate nodes
Parameter = establishModel(InitialParameter, 'gridFineness', 'high');   % denser intermediate nodes

% Manual mesh
manualGrid{1} = [1,2,1,7,1,1,3]; % shaft 1: sub-divisions between consecutive key nodes
manualGrid{2} = [1,3,4,3];       % shaft 2
Parameter = establishModel(InitialParameter, 'gridFineness', manualGrid);
```

<img src="readme/example2/MeshResultOfShaft1_manual_mesh.png" width="100%" height="" />
<img src="readme/example2/MeshResultOfShaft2_manual_mesh.png" width="100%" height="" />

#### Disable Auto-Plots

```matlab
Parameter = establishModel(InitialParameter, 'isPlotModel', false, 'isPlotMesh', false);
```

---

### 6.3 Time-Domain Analysis

Three steps after `establishModel`:

```matlab
generateDynamicEquation(Parameter);             % generate ODE function file
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14);  % integrate
signalProcessing(q, dq, t, Parameter, [25,35], 2^14, SwitchFigure); % post-process
```

Outputs:
- `q` — `[dofNum × nTimePoints]` displacement history
- `dq` — `[dofNum × nTimePoints]` velocity history
- `t` — `[1 × nTimePoints]` time vector

#### Solvers

| Name | Type | Recommended for |
|---|---|---|
| `'RK'` (default) | Fixed-step RK4 | Linear or mildly nonlinear systems |
| `'ode15s'` | Variable-step, implicit | **Hertzian contact, stiff systems** |
| `'ode45'` | Variable-step, explicit | Moderately nonlinear |
| `'ode23s'` | Variable-step, implicit | Alternative stiff solver |

```matlab
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'calculateMethod', 'ode15s');
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'calculateMethod', 'ode15s', ...
             'options', odeset('RelTol', 1e-6));
```

#### Initial Conditions

```matlab
% Zero (default)
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14);

% Custom displacement vector
q0 = zeros(Parameter.Mesh.dofNum, 1);
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, q0);

% Start from static equilibrium (auto-saved; recalculate with isFreshInitial)
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'isUseBalanceAsInitial', true);
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'isUseBalanceAsInitial', true, ...
             'isFreshInitial', true);
```

#### Downsampling

```matlab
downsampling = 5;
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'reduceInterval', downsampling);
% Pass the same factor to signalProcessing:
signalProcessing(q, dq, t, Parameter, [25,35], 2^14, SwitchFigure, 'reduceInterval', downsampling)
```

#### Suppress Status Diagram

```matlab
[q, dq, t] = calculateResponse(Parameter, [0,50], 2^14, 'isPlotStatus', false);
```

---

### 6.4 Frequency-Domain Analysis

Frequency-domain analysis works directly on `Parameter` after `establishModel` — **no code generation or time integration needed**.

#### Campbell Diagram

```matlab
max_rpm   = 10000;
exciteRad = linspace(1, max_rpm/60*2*pi, 500);  % rad/s sweep

[eigMatrix, criticalSpeeds] = calculateCampbell(Parameter, exciteRad, ...
    'isPlot',          true, ...
    'isFilter',        true, ...
    'filterMethod',    'slope', ...  % 'slope' or 'MAC'
    'isUseGyroMatrix', true);
```

| Output | Size | Description |
|---|---|---|
| `eigMatrix` | `[nModes × nSpeeds]` | Natural frequencies (rad/s) at each speed |
| `criticalSpeeds` | `[1 × K]` | Speeds where a mode line crosses the 1× excitation line |

| Option | Default | Description |
|---|---|---|
| `isPlot` | `false` | Draw the Campbell diagram |
| `isFilter` | `false` | Enable mode-tracking to prevent line crossings |
| `filterMethod` | `'MAC'` | Mode tracking: `'MAC'` (eigenvector correlation) or `'slope'` (frequency gradient) |
| `tolRad` | `0.1` | Rigid-body mode rejection threshold (rad/s) |
| `isUseGyroMatrix` | `true` | Include gyroscopic matrix in eigenvalue problem |

---

#### Mode Shapes

Visualize deformed rotor geometry at critical speeds:

```matlab
[ModeShapes, ZCoords] = calculateModeShape(Parameter, criticalSpeeds, ...
    'isPlot',          true, ...
    'direction',       'X', ...    % 'X' (horizontal) or 'Y' (vertical)
    'isUseGyroMatrix', true);
```

| Output | Description |
|---|---|
| `ModeShapes` | `{shaftNum × nSpeeds}` cell — displacement at each node per shaft |
| `ZCoords` | `{shaftNum × 1}` cell — global Z-axis positions of each shaft's nodes |

| Option | Default | Description |
|---|---|---|
| `isPlot` | `false` | Render deformed shaft geometry |
| `direction` | `'X'` | Extraction direction: `'X'` or `'Y'` |
| `scaleFactor` | `'auto'` | Auto scales max displacement to 1.5× outer shaft radius |
| `isUseGyroMatrix` | `true` | Include gyroscopic effects |

---

#### Unbalance Response

Solve the synchronous FRF across a speed sweep:

```matlab
% Single condition (from Parameter.Status)
response = calculateUnbalanceResponse(Parameter);
```

```matlab
% Multi-speed sweep — speedMatrix: [shaftNum × numConditions] (rad/s)
speedMatrix = [100, 200, 300;   % shaft 1 speed at each condition
               200, 400, 600];  % shaft 2 speed at each condition
response = calculateUnbalanceResponse(Parameter, speedMatrix);
```

```matlab
% Sweep with orbit trajectory plots
response = calculateUnbalanceResponse(Parameter, speedMatrix, ...
    'isPlot',           true, ...
    'plotNodeID',       [1, 4], ...     % nodes to visualize
    'plotConditionIdx', [50, 100]);     % speed-column indices to plot
```

Output `response` is a complex 3-D array `[dofNum × numConditions × shaftNum]`:

```
response(:, i, k)           → DOF amplitudes at condition i, excited by shaft k only
abs(response(dof, :, k))    → amplitude vs. speed (Bode magnitude) for DOF dof
angle(response(dof, :, k))  → phase vs. speed
```

| Option | Default | Description |
|---|---|---|
| `isPlot` | `false` | Plot steady-state orbit trajectories |
| `plotNodeID` | `1` | Node IDs to plot |
| `plotConditionIdx` | `1` | Speed conditions (columns of speedMatrix) to plot |

---

#### Speed-Dependent Bearings

When bearing stiffness/damping varies with speed (e.g., linearized Hertzian coefficients), use `inputSpeedDependentBearing*`:

```matlab
InitialParameter = inputEssentialParameterSingle2();
InitialParameter = inputSpeedDependentBearingSingle2(InitialParameter);

% Zero-out the standard Bearing struct (speed-dependent bearings replace it)
InitialParameter.Bearing.amount   = 0;
InitialParameter.Bearing.stiffness = [];
% ... zero out all Bearing fields ...

Parameter = establishModel(InitialParameter);
```

Inside `inputSpeedDependentBearingSingle2.m`, define the lookup table:

```matlab
SpeedDependentBearing.amount                  = 2;
SpeedDependentBearing.inShaftNo               = [1; 1];
SpeedDependentBearing.dofOfEachNodes          = [0; 0];          % must be 0
SpeedDependentBearing.positionOnShaftDistance = [0; 517.2]*1e-3; % m

SpeedDependentBearing.speed             = [100, 200, 300];        % rad/s lookup points

% Stiffness tables: [n bearings × m speed points] (N/m)
SpeedDependentBearing.stiffness         = [1e6, 5e6, 1e7; 2e6, 7e6, 2e7]; % Kxx
SpeedDependentBearing.stiffnessVertical = [...];  % Kyy
SpeedDependentBearing.stiffnessHV       = [...];  % Kxy (cross-coupled)
SpeedDependentBearing.stiffnessVH       = [...];  % Kyx (cross-coupled)

% Damping tables: [n bearings × m speed points] (Ns/m)
SpeedDependentBearing.damping           = [1e3, 2e3, 3e3; 4e3, 5e3, 6e3]; % Cxx
SpeedDependentBearing.dampingVertical   = [...];  % Cyy
SpeedDependentBearing.dampingHV         = [...];  % Cxy
SpeedDependentBearing.dampingVH         = [...];  % Cyx

SpeedDependentBearing.mass      = zeros(2,1); % must be 0
SpeedDependentBearing.isHertzian = zeros(2,1); % must be 0
```

`calculateUnbalanceResponse` automatically detects `ComponentSwitch.hasSpeedDependentBearing = true` and interpolates the lookup tables at each frequency step.

---

### 6.5 Post-Processing

```matlab
SwitchFigure.displacement     = true;  % time history (X, Y per node)
SwitchFigure.axisTrajectory   = true;  % 2-D orbit (X vs Y)
SwitchFigure.axisTrajectory3d = true;  % 3-D orbit (X vs Y vs t)
SwitchFigure.phase            = true;  % phase portrait (displacement vs velocity)
SwitchFigure.fftSteady        = true;  % steady-state FFT spectrum
SwitchFigure.fftTransient     = true;  % STFT waterfall / spectrogram
SwitchFigure.poincare         = true;  % Poincaré section
SwitchFigure.poincare_phase   = true;  % Poincaré + phase portrait overlay
SwitchFigure.saveFig          = true;  % save .fig (PNG always saved)
SwitchFigure.saveEps          = false; % save .eps vector graphics

% time_span: analysis window within the full simulation
signalProcessing(q, dq, t, Parameter, [25,35], 2^14, SwitchFigure)
```

All figures are saved to `signalProcess/`. For full name-value pair documentation, right-click `signalProcessing` in MATLAB and select **Help**.

If downsampling was used in `calculateResponse`, pass the same factor here:

```matlab
downsampling = 5;
signalProcessing(q, dq, t, Parameter, [25,35], 2^14, SwitchFigure, 'reduceInterval', downsampling)
```

---

## 7. References

- Haopeng Zhang, Steven Chatterton, Kaifu Zhang, Donglin Li, Runhan Li, Jin Chen, Andrea Riva, Shuai Gao, Kuan Lu and Paolo Pennacchi. *A Dynamic Simulation Tool for Multi-Spool Rotor Systems.* (Submitted)

- Haopeng Zhang, Runhan Li, Kuan Lu, Xiaohui Gu, Ruijuan Sang, Donglin Li. *Dynamic Behavior of Twin-Spool Rotor-Bearing System with Pedestal Looseness and Rub Impact.* Applied Sciences. 2024; 14(3):1181. [https://doi.org/10.3390/app14031181](https://doi.org/10.3390/app14031181)

---

*Copyright © 2021–2026 Haopeng Zhang, Northwestern Polytechnical University / Politecnico di Milano. Licensed under the [MIT License](LICENSE).*