%-----------------------------------------------------------------------------------
% This example creates a twin-spool rotor system with speed ratio 1.3 during run-up.
% This example shows new features to calculate unbalance response in the frequency domain.
%-----------------------------------------------------------------------------------
clc
clear
close all

% 1. Input all positional, physical, and geometric parameters of the rotor
InitialParameter = inputEssentialParameterTwinSpool(); % all shaft, disk, operation parameters are saved in this function file; you can open and check
InitialParameter = inputBearingHertzTwinSpool(InitialParameter); % all bearing parameters are saved in this file (But Hertzian force is not used in this example. You can control in this file)
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter); % all inter-shaft bearing parameters are saved in this file

% 2. Establish model
% Parameter = establishModel(InitialParameter); % this function will create the global matrices in workspace (check Parameter), model and mesh diagram in folders: <modelDiagram>, <meshDiagram>
% establish model manually
manualGrid{1} = [1,2,1,7,1,1,3]; % for shaft 1
manualGrid{2} = [1,3,4,3];       % for shaft 2
Parameter = establishModel(InitialParameter, "gridFineness", manualGrid, "isPlotMesh", false, "isPlotModel", false);


% =========================================================================
% Usage 1: Calculate the steady-state unbalance response using default settings
% =========================================================================
% By providing only 'Parameter', the function automatically calculates the 
% response for a SINGLE speed condition based on the preset struct data:
%   - Shaft 1 speed = Parameter.Status.vmax
%   - Shaft 2 speed = Parameter.Status.vmax * Parameter.Status.ratio
%
% Output 'response' is a Complex 3D array of size [dofNum, numSpeeds, shaftNum].
% In this case, size is [dofNum, 1, 2]:
%   - response(:, 1, 1): Vibration of the entire system excited EXCLUSIVELY by 
%                        the unbalance on Shaft 1 (at Shaft 1's speed frequency).
%   - response(:, 1, 2): Vibration of the entire system excited EXCLUSIVELY by 
%                        the unbalance on Shaft 2 (at Shaft 2's speed frequency).
% (Note: In the calculation, the total Gyroscopic matrix is accurately constructed 
%        by superimposing the gyroscopic effects of both spinning shafts.)

response = calculateUnbalanceResponse(Parameter);


% =========================================================================
% Usage 2: Calculate responses across a custom multi-condition speed sweep
% =========================================================================
% You can manually define a 'speedMatrix' of size [shaftNum, numSpeeds].
% Each COLUMN represents a distinct operating condition (speeds in rad/s).
speedMatrix = [100, 200, 300;   % Row 1: Shaft 1 speeds for Condition 1, 2, 3
               200, 400, 600];  % Row 2: Shaft 2 speeds for Condition 1, 2, 3

% Output 'response2' will be of size [dofNum, 3, 2].
% For instance, response2(:, 2, 1) stores the complex response at Condition 2 
% (where Shaft 1 is at 200 rad/s and Shaft 2 is at 400 rad/s), considering 
% ONLY the excitation forces from Shaft 1.
% This is extremely useful for extracting data to plot Bode or Waterfall diagrams.

response2 = calculateUnbalanceResponse(Parameter, speedMatrix);


% =========================================================================
% Usage 3: Calculate and automatically plot time-domain orbit trajectories
% =========================================================================
% Using modern Name-Value pair arguments, you can directly visualize the 
% steady-state orbit (Lissajous-like precessions) for specific nodes.
% 
% - "isPlot"           : true to enable trajectory plotting (default is false).
% - "plotNodeID"       : array of node IDs to visualize (e.g., [1, 2]).
% - "plotConditionIdx" : array of condition indices from speedMatrix to plot.
%                        (e.g., [2, 3] means plotting the 2nd and 3rd columns).
%
% This will generate 4 separate figures (Node 1 @ Cond 2, Node 1 @ Cond 3, 
% Node 2 @ Cond 2, Node 2 @ Cond 3). The plots automatically superimpose 
% the physical (real) time-domain displacements from all different frequency 
% excitation sources to form the true complex trajectory.

response3 = calculateUnbalanceResponse(Parameter, speedMatrix, ...
    "isPlot", true, ...
    "plotNodeID", [1, 2], ...
    "plotConditionIdx", [2, 3]);