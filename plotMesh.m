%% plotMesh - Visualize finite element mesh discretization for rotor systems
% Generates detailed schematic diagrams of mesh configurations with component annotations.
%
%% Syntax
%   plotMesh(Parameter)
%
%% Description
% |plotMesh| creates comprehensive visualizations of FEM meshes with:
% * Key component location markers
% * Hierarchical node classification
% * Automated legend generation
% * Publication-quality formatting
% Output includes:
% * Per-shaft diagram files (.fig/.png)
% * Standardized output directory ('meshDiagram')
%
%% Input Arguments
% * |Parameter| - System configuration structure containing:
%   * |Mesh|: Finite element mesh data [struct]
%     .keyPointsDistance   % Key node positions per shaft [nShafts×1 cell]
%     .nodeDistance        % Node positions per shaft [nShafts×1 cell]
%     .Node                % Node properties array [nNodes×1 struct]:
%       .name              % Node ID
%       .onShaftNo         % Parent shaft index
%       .onShaftDistance   % Axial position [m]
%       .diskNo            % Disk association ID
%       .bearingNo         % Bearing association ID
%       .isBearing         % Bearing node flag
%   * |Shaft|: Geometric parameters [struct array]
%     .amount              % Number of shafts
%     .totalLength         % Lengths [m] [nShafts×1]
%   * |ComponentSwitch|: Component activation status [struct]
%     .hasRubImpact        % Rub-impact status
%     .hasIntermediateBearing % Intermediate bearing status
%     .hasLoosingBearing   % Loose bearing status
%     .hasCouplingMisalignment % Coupling status
%     .hasCustom           % Custom component status
%
%% Output
% Generates in './meshDiagram/' directory:
%   MeshResultOfShaft[n].fig  % MATLAB figure file
%   MeshResultOfShaft[n].png  % Image file
% Where [n] = shaft index (1,2,...,N)
%
%% Visualization Features
% 1. Node Classification:
%    - Key nodes: Red circles
%    - Regular nodes: Black vertical ticks
%    - Bearing mass nodes: Purple markers
% 2. Component Markers:
%    - Disks: 'D#'
%    - Bearings: 'B#'
%    - Loose bearings: 'B# Loose'
%    - Intermediate bearings: 'IB#'
%    - Rub impacts: 'Rub#'
%    - Coupling misalignments: 'CpMis#'
% 3. Intelligent Layout:
%    - Automatic y-offset for overlapping bearings
%    - Context-aware legend generation
% 4. Formatting Standards:
%    - IEEE-compliant font (Times New Roman)
%    - Optimized dimensions for publications
%
%% Annotation System
% 1. Node Identification:
%    - Node IDs below shaft line
% 2. Component Tagging:
%    - Component tags above shaft line
% 3. Legend Categories:
%    * Geometric Markers:
%       - Key node ○
%       - Regular node │
%       - Bearing mass △
%       - Intermediate bearing ◊
%    * Component Tags:
%       - D: Disk
%       - B: Bearing
%       - IB: Inter-shaft bearing
%       - Rub: Rub impact
%       - CpMis: Coupling misalignment
%       - B Loose: Loose bearing
%
%% Example
% % Generate and view mesh diagrams
% sysConfig = meshModel(baseParams);
% plotMesh(sysConfig); 
% winopen('meshDiagram/MeshResultOfShaft1.png');
%
%% Implementation Details
% 1. Directory Management:
%   * Creates 'meshDiagram' directory
%   * Clears previous outputs
% 2. Per-Shaft Processing:
%   * Generates separate figures for each shaft
% 3. Layer Creation:
%   * Main shaft visualization layer
%   * Dedicated legend layer
% 4. Dynamic Positioning:
%   * Automatic spacing for multi-bearing nodes
%   * Adaptive vertical scaling
% 5. Output Optimization:
%   * Resolution-independent vector formats (.fig)
%   * Publication-ready raster formats (.png)
%
%% See Also
% meshModel, plot2DStandard, femShaft
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


function plotMesh(Parameter)
Shaft = Parameter.Shaft;
keyPointsDistance = Parameter.Mesh.keyPointsDistance;
nodeDistance = Parameter.Mesh.nodeDistance;
Node = Parameter.Mesh.Node;
ShaftElement = Parameter.Mesh.ShaftElement; 
%%
% Generate folder to save figures
hasFolder = exist('meshDiagram','dir');
if hasFolder
    delete meshDiagram/*.fig;
    delete meshDiagram/*.png;
else
    mkdir('meshDiagram');
end
%%
% Plot the key points on each shaft
for iShaft = 1:1:Shaft.amount
    
    % [Modified] Pre-calculate the maximum radius across BOTH modes for unified scaling.
    % This ensures that when switching between Geometric and Stiffness plots,
    % the visual scale remains constant, allowing users to observe radius differences.
    currentShaftElems = ShaftElement([ShaftElement.ShaftNo] == iShaft);
    
    % Extract all outer radii (Geometric)
    all_Ro_Geo = [currentShaftElems.outerRadius];
    % Extract all outer radii (Stiffness)
    all_Ro_Stiff = [currentShaftElems.outerRadiusStiff];
    
    % Find global max radius for this shaft
    maxR_Geo = max(all_Ro_Geo);
    if isempty(all_Ro_Stiff) || all(isnan(all_Ro_Stiff))
        maxR_Stiff = 0;
    else
        maxR_Stiff = max(all_Ro_Stiff);
    end
    
    % Determine the unified reference radius
    maxR_Unified = max([maxR_Geo, maxR_Stiff]);
    if maxR_Unified == 0, maxR_Unified = 0.05; end % Fallback
    
    % Calculate text gap based on the unified radius to keep labels stable
    text_gap = max(maxR_Unified * 0.2, 0.005); 

    % Define modes
    modes = {'Geometric', 'Stiffness'};
    
    for iMode = 1:length(modes)
        currentMode = modes{iMode};
        
        if strcmp(currentMode, 'Geometric')
            fileSuffix = '_Geometric';
            titleSuffix = ' (Geometric)';
            getRo = @(e) e.outerRadius;
            getRi = @(e) e.innerRadius;
            fillColor = [0.9, 0.9, 0.9]; % Light Gray
        else
            fileSuffix = '_Stiffness';
            titleSuffix = ' (Stiffness)';
            getRo = @(e) e.outerRadiusStiff;
            getRi = @(e) e.innerRadiusStiff;
            fillColor = [0.85, 0.9, 0.95]; % Light Blue-Gray
        end
        
        % Create figure
        figureName = ['Mesh Result of Shaft ', num2str(iShaft), titleSuffix];
        
        % [Modified] Removed 'Color', 'w'. Allow default background color.
        % Units set to normalized for responsive sizing.
        h = figure('name', figureName, 'Visible', 'off', ...
                   'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
        
        % Create Axes
        % Main Axes: Top 70%
        mainAx = axes('Parent', h, 'Units', 'normalized', 'Position', [0.05, 0.22, 0.9, 0.70]);
        % Legend Axes: Bottom 15%
        bottomAx = axes('Parent', h, 'Units', 'normalized', 'Position', [0.05, 0.02, 0.9, 0.15]);
        
        % [Modified] Enclose the main axes with a box (top and right borders)
        box(mainAx, 'on'); 
        
        % Plot shaft elements (Pipe walls)
        % Note: We iterate again for drawing, using specific mode data
        for k = 1:length(currentShaftElems)
            elem = currentShaftElems(k);
            L = elem.Length;
            Ro = getRo(elem);
            Ri = getRi(elem);
            
            % Find start index and position
            startNodeIdx = elem.NodeNo(1);
            zStart = Node(startNodeIdx).onShaftDistance;
            
            % Plot upper part (Ri to Ro)
            rectangle(mainAx, 'Position', [zStart, Ri, L, Ro-Ri], ...
                  'FaceColor', fillColor, 'EdgeColor', [0.2, 0.2, 0.2]); hold(mainAx, 'on');
            
            % Plot lower part (-Ro to -Ri)
            rectangle(mainAx, 'Position', [zStart, -Ro, L, Ro-Ri], ...
                  'FaceColor', fillColor, 'EdgeColor', [0.2, 0.2, 0.2]); hold(mainAx, 'on');
        end
        
        % -------------------------------------------------------------------
        
        % Plot centerline
        plotLine = plot(mainAx, [0 keyPointsDistance{iShaft}(end)] , [0 0]); hold(mainAx, 'on');
        plotLine.LineWidth = 1;      
        plotLine.Color = '#000000';  
        plotLine.LineStyle = '--';  
        
        % Plot key points
        plotKeyPoints = scatter(mainAx, keyPointsDistance{iShaft},...
                                zeros(length(keyPointsDistance{iShaft}), 1)); hold(mainAx, 'on');
        plotKeyPoints.SizeData = 45;
        plotKeyPoints.MarkerFaceColor = '#CA3636';
        plotKeyPoints.Marker = 'o';
        
        % Plot regular nodes
        nodeWithoutKeyPoints = setdiff(nodeDistance{iShaft}, keyPointsDistance{iShaft});
        plotNodes = scatter(mainAx, nodeWithoutKeyPoints,...
                            zeros(length(nodeWithoutKeyPoints), 1)); hold(mainAx, 'on');
        plotNodes.SizeData = 45;
        plotNodes.MarkerFaceColor = '#000000';
        plotNodes.MarkerEdgeColor = '#000000';
        plotNodes.Marker = '|';
        plotNodes.LineWidth = 1.5;
        
        % Find the nodes locating on iShaft
        nodeNum = length(Node);
        condition1 = zeros(1,nodeNum);
        for iNode = 1:1:nodeNum
            condition1(iNode) = ismember(iShaft,Node(iNode).onShaftNo);
        end 
        condition = condition1 & ([Node.isBearing] == false);
        NodeSegment = Node( condition ); % Nodes on iShaft
        segmentNum = length(NodeSegment);
        
        % Mark node IDs
        nodeName = cell(1,segmentNum);
        for iSegment = 1:1:segmentNum
                nodeName{iSegment} = num2str( NodeSegment(iSegment).name );
        end
        xText = [NodeSegment.onShaftDistance];
        
        % [Modified] Use maxR_Unified for positioning to keep labels stable across modes
        yText = -(maxR_Unified + text_gap) * ones(1,segmentNum); 
        text(mainAx, xText,yText,nodeName, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 8);
         
        % Mark element components
        elementName = cell(1,segmentNum);
        for iSegment = 1:1:segmentNum
            % Disk
            diskName = [];
            if ~isempty(NodeSegment(iSegment).diskNo)
                diskName = ['D ',num2str(NodeSegment(iSegment).diskNo)];
            end
            
            % Bearing and loosening bearing
            bearingName = [];
            if ~isempty(NodeSegment(iSegment).bearingNo)
                if Parameter.ComponentSwitch.hasLoosingBearing
                    isLoosing = NodeSegment(iSegment).isLoosingBearing;
                    if isLoosing
                        bearingName = ['B ',num2str(NodeSegment(iSegment).bearingNo),' Loose'];
                    else
                        bearingName = ['B ',num2str(NodeSegment(iSegment).bearingNo)];
                    end 
                else
                    bearingName = ['B ',num2str(NodeSegment(iSegment).bearingNo)];
                end          
            end 
              
            % Intermediate bearing
            interBearingName = [];
            if Parameter.ComponentSwitch.hasIntermediateBearing
                if ~isempty(NodeSegment(iSegment).interBearingNo)
                    interBearingName = ['IB ',num2str(NodeSegment(iSegment).interBearingNo)];
                end
            end
              
            % Rub-impact
            rubImpactName = [];
            if Parameter.ComponentSwitch.hasRubImpact
                if ~isempty(NodeSegment(iSegment).rubImpactNo)
                    rubImpactName = ['Rub ',num2str(NodeSegment(iSegment).rubImpactNo)];
                end
            end
             
            % Coupling misalignment
            couplingMisName = [];
            if Parameter.ComponentSwitch.hasCouplingMisalignment
                if ~isempty(NodeSegment(iSegment).couplingNo)
                    couplingMisName = ['CpMis ',num2str(NodeSegment(iSegment).couplingNo)];
                end
            end
            
            % Custom node
            customName = [];
            if Parameter.ComponentSwitch.hasCustom
                if ~isempty(NodeSegment(iSegment).customNo)
                    customName = ['Cus ',num2str(NodeSegment(iSegment).customNo)];
                end
            end
               
            % Assembling element names
            elementName{iSegment} = {   customName;...
                                        couplingMisName;...
                                        rubImpactName;...
                                        interBearingName;...
                                        bearingName;...
                                        diskName};
                                    
            % Delete empty values
            nullIndex = cellfun(@isempty, elementName{iSegment});
            notNullIndex = ~nullIndex;
            elementName{iSegment} = elementName{iSegment}(notNullIndex);
        end 
        
        xText = [NodeSegment.onShaftDistance];
        
        % [Modified] Use maxR_Unified for positioning
        yText = (maxR_Unified + text_gap) * ones(1,segmentNum);
        text(mainAx, xText,yText,elementName, 'HorizontalAlignment', 'center',...
             'VerticalAlignment', 'bottom', 'FontSize', 8, 'Interpreter', 'none');
        
        % Mark bearing nodes
        condition1 = zeros(1,nodeNum);
        for iNode = 1:1:nodeNum
            condition1(iNode) = ismember(iShaft,Node(iNode).onShaftNo);
        end 
        condition = condition1 & ([Node.isBearing] == true);
        NodeSegmentB = Node( condition ); % Nodes on iShaft
        segmentNumB = length(NodeSegmentB);
        nodeNameB = cell(1,segmentNumB); % Mark bearing node name
        posRecorder = zeros(segmentNumB,1);
        noInColumnRecorder = zeros(segmentNumB,1);
        
        for iSegment = 1:1:segmentNumB
            nodeNameB{iSegment} = num2str( NodeSegmentB(iSegment).name ); % Save node name
            % Judge type of bearing
            condition1 = Parameter.ComponentSwitch.hasIntermediateBearing; 
            condition2 = ~isempty(NodeSegmentB(iSegment).interBearingNo); 
            isInterBearing = condition1 && condition2;
            
            if isInterBearing
                % Find which shaft the InterBearing locates on
                shaftIndex = find(NodeSegmentB(iSegment).onShaftNo == iShaft); 
                xText = NodeSegmentB(iSegment).onShaftDistance(shaftIndex); 
            else
                xText = NodeSegmentB(iSegment).onShaftDistance;
            end
            
            % Calculate the y-position (stacking logic)
            if iSegment==1
                noInColumn = 1; 
            else
                previousNum = length( find(posRecorder(1:iSegment-1)==xText) );
                noInColumn = previousNum + 1;
            end
            posRecorder(iSegment) = xText; 
            noInColumnRecorder(iSegment) = noInColumn;
            
            % [Modified] Use maxR_Unified for positioning
            baseY = -(maxR_Unified + text_gap * 4); 
            yText = baseY - (noInColumn-1) * (text_gap * 3); 
            xMark = xText;
            yMark = yText;
            
            % Plot node name
            text(mainAx, xText,yText,nodeNameB{iSegment}, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', '#7E2F8E');
            
            % Plot marker
            plotBearingNodes = scatter(mainAx, xMark,yMark+text_gap, 'Parent', mainAx); hold(mainAx, 'on');
            plotBearingNodes.SizeData = 36;
            plotBearingNodes.MarkerFaceColor = '#7E2F8E';
            plotBearingNodes.MarkerEdgeColor = '#7E2F8E';
            if isInterBearing
                plotBearingNodes.Marker = 'd';
            else
                plotBearingNodes.Marker = '^';
            end
        end
         
        % Set axes limits and title
        axis(mainAx, 'equal'); 
        
        maxStackCount = max([1;noInColumnRecorder]);
        % [Modified] Calculate limits based on Unified Radius for consistent zooming
        upperLim = maxR_Unified + text_gap * 6; 
        lowerLim = -(maxR_Unified + text_gap * 4 + maxStackCount * (text_gap * 3)); 
        
        xlim(mainAx, [-0.05*Shaft.totalLength(iShaft), Shaft.totalLength(iShaft)*1.05]);
        ylim(mainAx, [lowerLim, upperLim]);
        
        title(mainAx, figureName);
        
        % Set legend axes position
        bottomAx.XTick = [];
        bottomAx.YTick = [];
        bottomAx.XTickLabel = [];
        bottomAx.YTickLabel = [];
        bottomAx.XColor = 'none';
        bottomAx.YColor = 'none';
        bottomAx.Color = 'none';
        bottomAx.XLim = [0, 1];
        bottomAx.YLim = [0, 1];
        hold(bottomAx, 'on');
        
        % Determine bearing existence for legend
        condition1 = [Node.isBearing];
        if isfield(Node, 'interBearingNo')
            n = numel(Node);
            values = cell(n, 1);
            for i = 1:n
                values{i} = Node(i).interBearingNo;
            end
            boolVector = double(~cellfun(@isempty, values));
            condition2 = boolVector(:)';
        else
            condition2 = boolean(zeros(size(condition1)));
        end
        hasOrdinaryBearing = sum(condition1 & ~condition2);
        hasInterBearing = sum(condition1 & condition2);
        
        % Assemble legend text
        textLegends = {};
        hasDisk = any(~cellfun(@isempty, {NodeSegment.diskNo}));
        hasBearing = any(~cellfun(@isempty, {NodeSegment.bearingNo}));
        hasLoosingBearing = Parameter.ComponentSwitch.hasLoosingBearing;
        hasInterBearingText = Parameter.ComponentSwitch.hasIntermediateBearing;
        hasRubImpact = Parameter.ComponentSwitch.hasRubImpact;
        hasCouplingMis = Parameter.ComponentSwitch.hasCouplingMisalignment;
        hasCustom = Parameter.ComponentSwitch.hasCustom;
        
        if hasDisk
            textLegends{end+1} = 'D: Disk';
        end
        if hasBearing
            textLegends{end+1} = 'B: Bearing';
        end
        if hasInterBearingText
            textLegends{end+1} = 'IB: Inter-shaft Bearing';
        end
        if hasCustom
            textLegends{end+1} = 'Cus: Customize Force';
        end
        if hasLoosingBearing
            textLegends{end+1} = 'B Loose: Loosening Bearing';
        end
        if hasRubImpact
            textLegends{end+1} = 'Rub: Rub Impact';
        end
        if hasCouplingMis
            textLegends{end+1} = 'CpMis: Coupling Misalignment';
        end
        
        % Integrate all text
        full_text = strjoin(textLegends, '       ');
        
        % Draw Legend Items
        xlim(bottomAx, [-1,12])
        xPositions = 0;
        x_pos_gap = 1.5;
        x_text_gap = 0.1;
        yPos = 0.75;
        
        % Key node legend
        s_key_point = scatter(bottomAx, xPositions, yPos);
        s_key_point.SizeData = 45;
        s_key_point.MarkerFaceColor = '#CA3636';
        s_key_point.Marker = 'o';
        s_key_point.MarkerEdgeColor = "none";
        text(bottomAx, xPositions+x_text_gap, yPos, 'Key Node', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
        xPositions = xPositions + x_pos_gap;
        
        % Node legend
        if ~isempty(nodeWithoutKeyPoints)
            s_node = scatter(bottomAx, xPositions, yPos);
            s_node.SizeData = 45;
            s_node.MarkerFaceColor = '#000000';
            s_node.MarkerEdgeColor = '#000000';
            s_node.Marker = '|';
            s_node.LineWidth = 1.5;
            text(bottomAx, xPositions+x_text_gap, yPos, 'Node', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
            xPositions = xPositions + x_pos_gap;
        end
        
        % Bearing legend
        if hasOrdinaryBearing
            s_beairng = scatter(bottomAx, xPositions,yPos);
            s_beairng.SizeData = 45;
            s_beairng.MarkerFaceColor = '#7E2F8E';
            s_beairng.MarkerEdgeColor = '#7E2F8E';
            s_beairng.Marker = '^';
            text(bottomAx, xPositions+x_text_gap, yPos, 'Mass at Bearing', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
            xPositions = xPositions + x_pos_gap + 0.5;
        end
        
        % Inter-bearing legend
        if hasInterBearing
            s_interbeairng = scatter(bottomAx, xPositions,yPos);
            s_interbeairng.SizeData = 45;
            s_interbeairng.MarkerFaceColor = '#7E2F8E';
            s_interbeairng.MarkerEdgeColor = '#7E2F8E';
            s_interbeairng.Marker = 'd';
            text(bottomAx, xPositions+x_text_gap, yPos, 'Mass at Inter-shaft Bearing', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
            xPositions = xPositions + x_pos_gap;
        end
        
        xPositionsText = 0;
        yPosText = 0.25; 
        
        text(bottomAx, xPositionsText, yPosText, full_text, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
            
        % Save figure
        set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
        figureName2 = ['meshDiagram/MeshResultOfShaft', num2str(iShaft), fileSuffix, '.fig'];
        savefig(h,figureName2)
        saveas(h, ['meshDiagram/MeshResultOfShaft', num2str(iShaft), fileSuffix, '.png'])
        close(h)
    
    end % end for iMode
end % end for iShaft
end % end function