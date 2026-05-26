# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**SpoolDyn** — an open-source MATLAB finite element rotor dynamics simulator for multispool systems (single or multi-shaft). It uses Timoshenko beam theory, lumped-parameter modeling, and code generation to simulate, analyze, and post-process rotor vibration responses.

- **Active development branch**: `develop` (16 commits ahead of `master`)
- **Published branch**: `master`
- **No build step, no package manager** — pure MATLAB `.m` files, run directly in MATLAB.

## Running Examples (Primary Way to Test)

Run any `main_example_N.m` script directly in MATLAB from the repository root. These serve as both usage demonstrations and integration tests:

| Script | What it shows |
|--------|--------------|
| `main_example_1.m` | Single-spool run-up, time-domain response |
| `main_example_2.m` | Twin-spool with Hertzian contact |
| `main_example_3.m` | Variable cross-section shaft, cross-term bearings |
| `main_example_4.m` | Twin-spool with intermediate (inter-shaft) bearing |
| `main_example_5_campbell_and_mode.m` | Campbell diagram + mode shapes |
| `main_example_6_solve_in_frequency_domain.m` | FRF / unbalance response sweep |
| `main_example_7_speed_independent_bearing.m` | Speed-dependent bearing in FRF solver |

To run a single example via MCP:
```matlab
run('main_example_1.m')
```

Output folders auto-created in the working directory: `modelDiagram/`, `meshDiagram/`, `signalProcess/`.

## Simulation Pipeline (Procedural, Always in This Order)

```
1. InitialParameter = input*()           % Define geometry, material, speed profile
2. Parameter = establishModel(InitialParameter)  % FEM assembly → Parameter.Matrix, Parameter.Mesh
3. generateDynamicEquation(Parameter)    % Code-generation → writes dynamicEquation.m to disk
4. [q, dq, t] = calculateResponse(...)  % Time-domain ODE integration
5. signalProcessing(q, dq, t, ...)      % Post-processing & plots
```

For **frequency-domain / modal** workflows, steps 3–5 are replaced:
```matlab
[eigMatrix, criticalSpeed] = calculateCampbell(Parameter, exciteRad)
[ModeShapes, ZCoords]      = calculateModeShape(Parameter, criticalSpeed)
response = calculateUnbalanceResponse(Parameter, speedMatrix)
```

## Core Data Structure: `Parameter` Struct

The single nested struct threaded through all functions:

| Field | Added by | Contents |
|-------|----------|----------|
| `Shaft`, `Disk`, `Bearing`, `IntermediateBearing` | `input*()` | Geometry, material, position |
| `Status` | `input*()` | Speed profile: `vmax`, `acceleration`, `ratio`, `isDeceleration`, `isUseCustomize` |
| `ComponentSwitch` | `input*()` | Feature flags: `hasGravity`, `hasHertzianForce`, `hasRubImpact`, `hasLoosingBearing`, `hasCouplingMisalignment`, `hasIntermediateBearing`, `hasCustom` |
| `Mesh` | `meshModel()` inside `establishModel` | Nodes, DOF intervals, element properties |
| `Matrix` | `establishModel()` | Sparse global matrices: `mass`, `stiffness`, `damping`, `gyroscopic`, `matrixN`, `gravity`, `unbalance`, `HerzianParameter` |

### Input Format Versions (V1 → V2)

Shaft properties use **cell arrays per segment** in V2 (current). V1 used scalar/vector double arrays. `standardizeInputParameter()` auto-converts V1 → V2 at the start of `establishModel`. When creating new `input*()` files, use V2 cell array format:
```matlab
Shaft.outerRadius = {[0.05; 0.06]; [0.04]};  % cell per shaft, vector per segment
```

## Key Architectural Decisions

### Code Generation
`generateDynamicEquation(Parameter)` writes `dynamicEquation.m` and `dynamicEquationNoMass.m` to the repository root at runtime. These files are **overwritten every run** and tailored to the specific system (hardcoded DOF indices, force terms). The `.gitignore` excludes them. Never edit these generated files manually.

### FEM Assembly
- Each shaft element: 8 DOF (4 per node: x, θx, y, θy)  
- Global matrices assembled via `addElementIn(A, B, position)` for submatrix insertion  
- Rayleigh damping: `C = α·M + β·K` applied to rotating parts only  
- All final matrices are `sparse`  
- Small values below `matrix_value_tol` (default `1e-12`) are zeroed

### Hertzian Contact
- Bearing Hertzian parameters collected by `collectHerzianParameter()` into `Matrix.HerzianParameter`  
- Force evaluated at each timestep by `hertzianForce()` (called inside generated `dynamicEquation.m`)  
- Jacobians: `hertzianForceJacobian()`, `hertzianForceLocalJacobian()` — enable `isUseJacobian=true` in `calculateResponse`

### Speed-Dependent Bearings (Develop branch)
`inputSpeedDependentBearingSingle2.m` defines bearing matrices as functions of speed. `updateSpeedDependentMatrices()` rebuilds matrices at each frequency step in the FRF solver.

### Multi-Shaft Speed Ratios
Shaft 1 is always the reference. Other shafts follow: `ω_n = vmax × ratio(n-1)`. The `Status.ratio` vector length = `Shaft.amount - 1`.

## Adding a New Configuration

1. Copy an existing `input*()` file as a template matching your topology (single/twin-spool, with/without inter-shaft bearing)
2. Set `ComponentSwitch` flags for the features you need
3. For Hertzian bearings: also call `inputBearingHertz*(InitialParameter)` to append bearing params
4. For inter-shaft bearings: also call `inputIntermediateBearing*(InitialParameter)`

## Important Utility Functions

- `addElementIn(A, B, pos)` — insert submatrix B into A at [row,col] position (used everywhere in FEM assembly)
- `assembleLinear(...)` — global matrix assembly  
- `standardizeInputParameter(p)` — V1→V2 format conversion (always called first in `establishModel`)
- `getDesignPalette()` — returns consistent color palette for all plots  
- `cell2string(c)` — converts cell arrays to strings (used inside code generation)
- `get_fft_components(signal, fs)` — extract FFT amplitude/phase at specific frequencies

## Output Folders

| Folder | Created by | Contents |
|--------|-----------|----------|
| `modelDiagram/` | `plotModel()` | Schematic diagrams |
| `meshDiagram/` | `plotMesh()` | FEM mesh visualizations |
| `signalProcess/` | `signalProcessing()` | Time histories, FFT, trajectories |

## Naming Conventions

- `input*()` functions: define system parameters, return `InitialParameter` or append to it  
- `fem*()` functions: generate component-level FEM matrices  
- `calculate*()` functions: perform analysis computations  
- `generate*()` functions: code generation or force generation  
- `plot*()` functions: visualization  
- Files prefixed `main_example_` are runnable scripts (not functions)
