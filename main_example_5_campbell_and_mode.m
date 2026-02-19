%-----------------------------------------------------------------------------------
% This example creates a twin-spool rotor system with speed ratio 1.3 during run-up.
% This example shows new features to get Campbell diagram and modes
%-----------------------------------------------------------------------------------

clc
clear
close all

% input all positional, physical, and geometric parameters of the rotor
InitialParameter = inputEssentialParameterTwinSpool(); % all shaft, disk, operation parameters are saved in this function file; you can open and check
InitialParameter = inputBearingHertzTwinSpool(InitialParameter); % all bearing parameters are saved in this file (But Hertzian force is not used in this example. You can control in this file)
InitialParameter = inputIntermediateBearingTwinSpool(InitialParameter); % all inter-shft bearing parameters are saved in this file

% establish model automatically
% Parameter = establishModel(InitialParameter); % this function will create the global matrices in workspace (check Parameter), model and mesh diagram in folders: <modelDiagram>, <meshDiagram>
% establish model manually
manualGrid{1} = [1,2,1,7,1,1,3]; % for shaft 1
manualGrid{2} = [1,3,4,3]; % for shaft 2
Parameter = establishModel(InitialParameter, "gridFineness", manualGrid);


% calculate campell diagram
max_rpm = 10000; % rpm
exciteRad = linspace(1,max_rpm/60*2*pi, 500);
[eigMatrix, criticalSpeed] = calculateCampbell(Parameter,exciteRad,"isPlot",true,"isFilter",true,"filterMethod","slope");

% calculate mode shape
[ModeShapes, ZCoords] = calculateModeShape(Parameter, criticalSpeed,"isPlot", true, "direction", "X");