%--------------------------------------------------------------------------
% This example creates a two-disk rotor system with time independent
% bearings (Just for solving in frequency domain)
%--------------------------------------------------------------------------

clc
clear
close all

% input all positional, physical, and geometric parameters of the rotor
InitialParameter = inputEssentialParameterSingle2(); % all input parameters are saved in this function file; you can open and check

% input parameters of speed independent bearings
InitialParameter = inputSpeedDependentBearingSingle2(InitialParameter);

% set the normal bearing parameters as zero
InitialParameter.Bearing.amount          = 0;
InitialParameter.Bearing.inShaftNo       = [];
InitialParameter.Bearing.dofOfEachNodes  = []; % if mass=0, dof must be 0 
InitialParameter.Bearing.positionOnShaftDistance = []; % m
InitialParameter.Bearing.stiffness       = []; % N*m
InitialParameter.Bearing.stiffnessVertical = []; % N*m
InitialParameter.Bearing.damping         = []; % N*s/m
InitialParameter.Bearing.dampingVertical = []; % N*s/m
InitialParameter.Bearing.mass            = []; % kg
InitialParameter.Bearing.isHertzian      = [];

% establish model automatically
Parameter = establishModel(InitialParameter); % this function will create the global matrices in workspace (check Parameter), model and mesh diagram in folders: <modelDiagram>, <meshDiagram>

% Calculate in frequency domain
speedMatrix = linspace(0,500, 100);
response3 = calculateUnbalanceResponse(Parameter, speedMatrix, ...
    "isPlot", true, ...
    "plotNodeID", [1, 4], ...
    "plotConditionIdx", [50, 100]);


%-----------------------------------------------------------------------------------------------------------
% After running this script, please check folders in the root: <meshDiagram> <modelDiagram> <signalProcess>
%-----------------------------------------------------------------------------------------------------------

