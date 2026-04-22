%% plotModel - Visualize 3D geometry of multi-shaft rotor systems
%
% This function generates detailed 3D schematic diagrams of rotor systems,
% including shafts, disks, bearings, and intermediate bearings. It visualizes
% both the physical geometric model and the equivalent stiffness model.
%
%% Syntax
%   plotModel(InitialParameter)
%
%% Description
% |plotModel| creates comprehensive 3D visualizations of rotor systems:
% * Renders cylindrical shafts with inner/outer diameters (segmented)
% * Displays disk geometries at specified positions
% * Visualizes bearing housings as triangular blocks
% * Automatically calculates intermediate bearing alignments
% * Generates diagrams for:
%   1. Geometric Model (Physical dimensions)
%   2. Stiffness Model (Equivalent stiffness diameters, if available)
%
%% Input Arguments
% * |InitialParameter| - System configuration structure containing:
%   * |Shaft|: [1×1 struct]              % Shaft properties
%     .amount             % Number of shafts [scalar]
%     .segmentLength      % Segment lengths [Cell array of vectors]
%     .outerRadius        % Physical outer radii [Cell array of vectors]
%     .innerRadius        % Physical inner radii [Cell array of vectors]
%     .outerRadiusStiff   % (Optional) Stiffness outer radii [Cell array]
%     .innerRadiusStiff   % (Optional) Stiffness inner radii [Cell array]
%   * |Disk|: [1×1 struct]               % Disk parameters
%     .amount             % Number of disks [scalar]
%     .inShaftNo          % Parent shaft indices [M×1 vector]
%     .positionOnShaftDistance % Axial positions [m] [M×1 vector]
%     .outerRadius        % Disk radii [m] [M×1 vector]
%     .thickness          % Disk thicknesses [m] [M×1 vector]
%   * |Bearing|: [1×1 struct]            % Bearing parameters
%     .amount             % Number of bearings [scalar]
%     .inShaftNo          % Parent shaft indices [K×1 vector]
%     .positionOnShaftDistance % Axial positions [m] [K×1 vector]
%   * |IntermediateBearing|: [1×1 struct] % (Optional) Intermediate bearings
%     .amount             % Number of intermediate bearings [scalar]
%     .betweenShaftNo     % Connected shaft pairs [L×2 matrix]
%     .positionOnShaftDistance % Connection positions [m] [L×2 matrix]
%
%% Output
% Creates in './modelDiagram' directory:
% * |diagramOfShaft[n].png|           % Physical shaft images
% * |theWholeModel.png|               % Composite physical system image
% * |diagramOfShaft[n]_Stiffness.png| % Stiffness model shaft images (if applicable)
% * |theWholeModel_Stiffness.png|     % Composite stiffness system image (if applicable)
% * corresponding |.fig| files for all outputs
%
%% Visualization Features
% 1. Dual Model Rendering:
%    * Geometric Mode: Visualizes actual physical dimensions
%    * Stiffness Mode: Visualizes equivalent stiffness diameters (useful for stepped shafts)
% 2. Component Rendering:
%    * Shafts: Segmented hollow cylinders
%    * Disks: Solid cylinders with specified thickness
%    * Bearings: Triangular housing structures
% 3. Automatic Alignment:
%    * Calculates shaft position offsets from intermediate bearings
%    * Maintains geometric relationships between connected shafts
% 4. Lighting and Rendering:
%    * Dual directional lighting (cool/warm tones)
%    * Gouraud shading for smooth surfaces
%
%% Implementation Details
% 1. Model Iteration:
%    * Loops through 'Geometric' and 'Stiffness' modes
%    * Selects appropriate radius data for each pass
% 2. Offset Calculation:
%    * Computes shaft position adjustments based on intermediate bearings
% 3. Directory Management:
%    * Creates 'modelDiagram' directory if missing
%    * Clears previous outputs before generation
% 4. Composite Diagram:
%    * Combines all shafts into single visualization
%    * Preserves component colors and styles
%
%% Example
% % Configure and visualize rotor system
% rotorParams = inputEssentialParameter(); % Load system parameters
% plotModel(rotorParams);                   % Generate diagrams
% % View composite diagram
% winopen('modelDiagram/theWholeModel.png');
%
%% Component Rendering Notes
% * Shaft Dimensions:
%    - Geometric Mode: Uses Shaft.outerRadius / Shaft.innerRadius
%    - Stiffness Mode: Uses Shaft.outerRadiusStiff / innerRadiusStiff
%    - If Stiffness radii are NaN/missing, falls back to geometric radii
% * Bearing Housing:
%    - Housing hole size adapts to the shaft diameter of the active mode
%
%% See Also
% addCylinder, addTriangularBlock, CombFigs, inputEssentialParameterBO
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


function plotModel(InitialParameter)
Shaft = InitialParameter.Shaft;
Disk = InitialParameter.Disk;
Bearing = InitialParameter.Bearing;

%% 1. Calculate the offset (Shared for both models)
% Offsets are based on segment lengths, which are assumed constant between 
% geometric and stiffness models.
offsetPosition = zeros(Shaft.amount,1);
if isfield(InitialParameter,'IntermediateBearing')
    InterBearing = InitialParameter.IntermediateBearing; 
    for iShaft = 1:1:Shaft.amount
        if iShaft == 1
            offsetPosition(1) = 0; 
        else
            for iInterBearing = 1:1:InterBearing.amount
                if InterBearing.betweenShaftNo(iInterBearing,2) == iShaft
                   basicShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 1);
                   laterShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 2);
                   basicShaftNo = InterBearing.betweenShaftNo(iInterBearing,1);
                   offsetPosition(iShaft) = basicShaftD - laterShaftD + offsetPosition(basicShaftNo);
                end 
            end 
        end 
    end 
end 

%% 2. Setup Output Directory
hasFolder = exist('modelDiagram','dir');
if hasFolder
    delete modelDiagram/*.fig;
    delete modelDiagram/*.png;
else
    mkdir('modelDiagram');
end

%% 3. Visualization Main Loop (Runs twice: Geometry & Stiffness)
% mode 1: Geometric Model (Physical)
% mode 2: Stiffness Model (Mathematical)
modelModes = {'Geometric', 'Stiffness'};
for iMode = 1:length(modelModes)
    currentMode = modelModes{iMode};
    
    % --- DATA SELECTION SWITCH ---
    if strcmp(currentMode, 'Geometric')
        % Use physical dimensions
        activeOuterRadius = Shaft.outerRadius;
        activeInnerRadius = Shaft.innerRadius;
        fileSuffix = ''; % Standard filename
        figTitlePrefix = 'Geometric Model';
    else
        % Use stiffness dimensions
        % Check if stiffness data exists, if not, skip or fallback
        if isfield(Shaft, 'outerRadiusStiff') && ~isempty(Shaft.outerRadiusStiff)
            activeOuterRadius = Shaft.outerRadiusStiff;
            activeInnerRadius = Shaft.innerRadiusStiff;
            fileSuffix = '_Stiffness';
            figTitlePrefix = 'Stiffness Model';
        else
            warning('Stiffness radius data not found. Skipping Stiffness Model plot.');
            continue; 
        end
    end
    
    % Create persistent figure for final composite of CURRENT mode
    wholeFig = figure('Visible', 'off', 'Name', [figTitlePrefix, ' - Whole System']);
    wholeAxes = axes(wholeFig);
    hold(wholeAxes, 'on');
    
    figureName = cell(Shaft.amount,1);
    for iShaft = 1:1:Shaft.amount
        h = figure('visible','off', 'Name', [figTitlePrefix, ' - Shaft ', num2str(iShaft)]);
        ax_h = axes(h);
        hold(ax_h, 'on'); 
        
        % --- DRAW SHAFT (Segmented) ---
        segLengths = Shaft.segmentLength{iShaft}; 
        % Select the active radius for the current mode
        segOuterR  = activeOuterRadius{iShaft}; 
        segInnerR  = activeInnerRadius{iShaft};
        
        currentLocalZ = 0; 
        
        for iSeg = 1:length(segLengths)
            L_seg = segLengths(iSeg);
            R_out = segOuterR(iSeg);
            R_in  = segInnerR(iSeg);
            
            % If Stiffness Radius is NaN (input convention), fallback to Geo Radius
            if isnan(R_out), R_out = Shaft.outerRadius{iShaft}(iSeg); end
            if isnan(R_in),  R_in  = Shaft.innerRadius{iShaft}(iSeg); end
            
            positionX = offsetPosition(iShaft) + currentLocalZ + L_seg/2;
            position = [positionX, 0, 0]; 
            
            NODES = 20;
            axisName = 'x';
            addCylinder(ax_h, position, R_out, R_in, L_seg, NODES, axisName);
            
            currentLocalZ = currentLocalZ + L_seg;
        end
        
        % --- DRAW DISK (Invariant) ---
        % Disks are physical masses, usually shown same in both models for reference
        for iDisk = 1:1:Disk.amount
            if Disk.inShaftNo(iDisk) == iShaft
                positionX = Disk.positionOnShaftDistance(iDisk) + offsetPosition(iShaft);
                position = [positionX, 0, 0];
                % Disks always use physical dimensions
                addCylinder(ax_h, position, Disk.outerRadius(iDisk), Disk.innerRadius(iDisk), ...
                            Disk.thickness(iDisk), 30, 'x');
            end 
        end 
        
        % --- DRAW BEARING (Standard) ---
        for iBearing = 1:1:Bearing.amount
            if Bearing.inShaftNo(iBearing) == iShaft
                positionX = Bearing.positionOnShaftDistance(iBearing) + offsetPosition(iShaft);
                position = [positionX, 0, 0];
                
                % 1. Find the shaft radius AT THIS MODE for the bearing hole
                bearingLoc = Bearing.positionOnShaftDistance(iBearing);
                localShaftR = 0;
                tempCursor = 0;
                for k = 1:length(segLengths)
                    if bearingLoc >= tempCursor && bearingLoc <= (tempCursor + segLengths(k))
                        localShaftR = segOuterR(k); % Use ACTIVE radius
                        if isnan(localShaftR), localShaftR = Shaft.outerRadius{iShaft}(k); end
                        break;
                    end
                    tempCursor = tempCursor + segLengths(k);
                end
                if localShaftR == 0, localShaftR = max(segOuterR); end 
                
                % 2. Housing size (Visual only)
                maxShaftR_Phys = max(Shaft.outerRadius{iShaft}); 
                if isempty(Disk.outerRadius)
                    maxDiskR = 0;
                else
                    maxDiskR = max(Disk.outerRadius); 
                end
                
                height = max(maxDiskR * 1.25, maxShaftR_Phys * 2.5);
                width = height;
                thickness = 0.01; 
                if ~isempty(Disk.thickness), thickness = min(Disk.thickness) * 0.6; end
                
                RotateInfo.isRotate = true;
                RotateInfo.oringin = [0,0,0];
                RotateInfo.direction = [1,0,0];
                RotateInfo.angle = 90;
                
                addTriangularBlock(ax_h, position, localShaftR, height, width, thickness, 15, 'x', RotateInfo);     
            end 
        end 

        % --- DRAW SPEED-DEPENDENT BEARING ---
        if isfield(InitialParameter, 'SpeedDependentBearing')
            SpdBearing = InitialParameter.SpeedDependentBearing;
            for iSpdBearing = 1:1:SpdBearing.amount
                if SpdBearing.inShaftNo(iSpdBearing) == iShaft
                    positionX = SpdBearing.positionOnShaftDistance(iSpdBearing) + offsetPosition(iShaft);
                    position = [positionX, 0, 0];
                    
                    % 1. Find the local shaft radius AT THIS MODE for the bearing hole
                    % This visualizes if the bearing is attached to the physical or stiffness radius
                    bearingLoc = SpdBearing.positionOnShaftDistance(iSpdBearing);
                    localShaftR = 0;
                    tempCursor = 0;
                    for k = 1:length(segLengths)
                        if bearingLoc >= tempCursor && bearingLoc <= (tempCursor + segLengths(k))
                            localShaftR = segOuterR(k); % Use ACTIVE radius
                            if isnan(localShaftR), localShaftR = Shaft.outerRadius{iShaft}(k); end
                            break;
                        end
                        tempCursor = tempCursor + segLengths(k);
                    end
                    if localShaftR == 0, localShaftR = max(segOuterR); end 
                    
                    % 2. Calculate visualization parameters (Housing dimensions)
                    maxShaftR_Phys = max(Shaft.outerRadius{iShaft}); 
                    maxDiskR = 0;
                    if ~isempty(Disk.outerRadius), maxDiskR = max(Disk.outerRadius); end
                    
                    height = max(maxDiskR * 1.25, maxShaftR_Phys * 2.5);
                    width = height;
                    thickness = 0.01; 
                    if ~isempty(Disk.thickness), thickness = min(Disk.thickness) * 0.6; end
                    
                    RotateInfo.isRotate = true;
                    RotateInfo.oringin = [0,0,0];
                    RotateInfo.direction = [1,0,0];
                    RotateInfo.angle = 90;
                    
                    % 3. Call the underlying drawing function
                    addTriangularBlock(ax_h, position, localShaftR, height, width, thickness, 15, 'x', RotateInfo);     
                end 
            end 
        end
        
        % --- Lighting and Saving Individual Shafts ---
        % Copy objects to composite figure
        allObjs = findall(h, 'type','axes');
        toCopy = allchild(allObjs);
        copyobj(toCopy, wholeAxes);
        
        light(ax_h, 'Position', [-1 -1 1], 'Color', [0.8 0.8 1]); 
        light(ax_h, 'Position', [1 1 1], 'Color', [1 0.9 0.8]);
        lighting(ax_h, 'gouraud');
        axis(ax_h, 'equal'); grid(ax_h, 'on'); view(ax_h, 3);
        title(ax_h, [figTitlePrefix, ': Shaft ', num2str(iShaft)]);
        
        % Save individual shaft figure with Suffix
        set(h,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
        figureName{iShaft} = ['modelDiagram/diagramOfShaft',num2str(iShaft), fileSuffix, '.fig'];
        savefig(h,figureName{iShaft})
        pngName = ['modelDiagram/diagramOfShaft',num2str(iShaft), fileSuffix, '.png'];
        saveas(h, pngName)
        close(h)
        
    end % End Shaft Loop
    
    % --- Finalize Composite Figure for this Mode ---
    view(wholeAxes, 3); 
    grid(wholeAxes, 'on');
    axis(wholeAxes, 'equal');
    xlabel(wholeAxes, 'X [m]'); ylabel(wholeAxes, 'Y [m]'); zlabel(wholeAxes, 'Z [m]');
    title(wholeAxes, [figTitlePrefix, ': Whole System']);
    light(wholeAxes, 'Position', [-1 -1 1], 'Color', [0.8 0.8 1]);
    light(wholeAxes, 'Position', [1 1 1], 'Color', [1 0.9 0.8]);
    lighting(wholeAxes, 'gouraud');
    
    % Save composite figure with Suffix
    set(wholeFig,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    savefig(wholeFig, ['modelDiagram/theWholeModel', fileSuffix, '.fig']);
    saveas(wholeFig, ['modelDiagram/theWholeModel', fileSuffix, '.png']);
    close(wholeFig);
end % End Mode Loop
end