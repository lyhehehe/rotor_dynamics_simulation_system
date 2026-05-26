%% meshModel - Generate discretized finite element mesh for multi-shaft rotor systems
%
% This function creates a comprehensive nodal mesh structure for finite element
% analysis of rotor-bearing systems. It handles automatic mesh generation based
% on shaft segmentation, identifies key components, assigns physical properties 
% to elements, and manages DOF mapping.
%
%% Syntax
%   Parameter = meshModel(InitialParameter)
%   Parameter = meshModel(InitialParameter, gridFineness)
%   Parameter = meshModel(InitialParameter, manualGrid)
%
%% Description
% |meshModel| performs sophisticated mesh generation for rotor dynamics models:
% * Automatically detects critical points (segment boundaries, disks, bearings)
% * Maps finite elements to specific shaft segments for property assignment
% * Implements adaptive mesh refinement strategies
% * Establishes node-component relationships
% * Manages degree-of-freedom (DOF) allocation
%
%% Input Arguments
% * |InitialParameter| - System configuration structure containing:
%   * |Shaft|: Shaft properties [struct array]
%     .totalLength         % Total length of each shaft [m]
%     .segmentLength       % Lengths of distinct segments [cell array]
%     .outerRadius         % Outer radii per segment [cell array]
%     .innerRadius         % Inner radii per segment [cell array]
%     .density             % Density per segment [cell array]
%     .eccentricity        % Eccentricity per segment [cell array]
%     .dofOfEachNodes      % DOFs per node
%   * |Disk|: Disk parameters [struct array]
%     .inShaftNo           % Parent shaft index
%     .positionOnShaftDistance % Axial position [m]
%   * |Bearing|: Bearing parameters [struct array]
%     .inShaftNo            % Parent shaft index
%     .positionOnShaftDistance % Axial position [m]
%   * |ComponentSwitch|: Component activation flags [struct]
%     .hasRubImpact, .hasIntermediateBearing, .hasLoosingBearing, etc.
%
% * |gridFineness| - Mesh resolution specification [string]:
%   * |'low'|: Coarse mesh (default, 1 element between key nodes)
%   * |'middle'|: Medium mesh (4 elements between key nodes)
%   * |'high'|: Fine mesh (10 elements between key nodes)
%
% * |manualGrid| - Custom mesh definition [cell array]:
%   {n} = [e1 e2 ... em]  % Element counts per segment for shaft n
%
%% Output Structure
% * |Parameter| - Enhanced system configuration with mesh data:
%   * |Mesh|: Discretization results [struct]
%     .nodeDistance        % Node positions per shaft [1×nShafts cell]
%     .Node                % Node properties [nNodes×1 struct]:
%       .name, .onShaftNo, .onShaftDistance, .dof, .isBearing, etc.
%     .ShaftElement        % Element-specific properties [nElements×1 struct]:
%       .NodeNo            % Global node indices [NodeL, NodeR]
%       .ShaftNo           % Parent shaft index
%       .Length            % Element length
%       .outerRadius       % Element outer radius
%       .innerRadius       % Element inner radius
%       .eccentricity      % Element eccentricity
%       .eccentricityPhase % Element eccentricity phase
%       .density, .elasticModulus, .poissonRatio
%     .dofInterval        % DOF ranges [nNodes×2]
%     .dofOnNodeNo         % Node mapping for DOFs [nDOFs×1]
%     .nodeNum             % Total node count
%     .dofNum              % Total DOF count
%   * Updated component fields with node mappings
%
%% Key Algorithms
% 1. Key Node Identification:
%    * Collects critical positions: Shaft segment boundaries, disks, bearings.
%    * Implements proximity merging (L/5000 threshold).
% 2. Mesh Segmentation:
%    * Subdivides shaft segments based on |gridFineness| or |manualGrid|.
% 3. Property Mapping:
%    * Matches generated finite elements to input shaft segments.
%    * Assigns specific geometric (radii) and material (density, E) 
%      properties to each element.
% 4. Node-Component Mapping:
%    * Associates physical components with nearest node.
% 5. Special Component Handling:
%    * Generates additional nodes for bearing masses.
%
%% Example
% % Medium-resolution automatic mesh
% sysConfig = meshModel(baseParams, 'middle');
%
% % Custom mesh for dual-shaft system
% manualGrid = {[3 2 4], [5 1]}; % Shaft1: 3|2|4 elements, Shaft2: 5|1
% sysConfig = meshModel(baseParams, manualGrid);
%
%% Implementation Notes
% * Node Merging Threshold: L/5000 (L = shaft length)
% * Element Properties: Each finite element inherits properties (including
%   eccentricity) from the specific shaft segment it resides in.
% * Bearing Nodes: Bearings with mass generate additional system nodes.
%
%% See Also
% establishModel
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


function Parameter = meshModel(varargin)   
% define the constant: this function will merge key nodes whose distance
% from other nodes to itself < totalLength/THRESHOLD_COEFFICIENT
THRESHOLD_COEFFICIENT = 5000; 

% --- Input Parsing ---
if isempty(varargin)
    error('Insufficient input parameters')
else
    InitialParameter = varargin{1};
end

switch length(varargin)
    case 1
        isAutoMesh = true;
        gridFineness = 'low';
    case 2
        if ischar(varargin{2})
            gridFineness = varargin{2};
            isAutoMesh = true;
        elseif iscell(varargin{2})
            manualGrid = varargin{2};
            isAutoMesh = false;
        else
            error('Please input char or cell data for the second parameter')
        end
    otherwise 
        error('too much input parameter')
end

Shaft = InitialParameter.Shaft;
Disk = InitialParameter.Disk;
Bearing = InitialParameter.Bearing;

% Component flags
if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasRubImpact') && InitialParameter.ComponentSwitch.hasRubImpact
    RubImpact = InitialParameter.RubImpact; hasRub = true;
else
    hasRub = false; 
end

if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasIntermediateBearing') && InitialParameter.ComponentSwitch.hasIntermediateBearing
    InterBearing = InitialParameter.IntermediateBearing; hasInterBearing = true;
else
    hasInterBearing = false; 
end

if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasCouplingMisalignment') && InitialParameter.ComponentSwitch.hasCouplingMisalignment
    Coupling = InitialParameter.CouplingMisalignment; hasCoupling = true;
else
    hasCoupling = false; 
end

if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasLoosingBearing') && InitialParameter.ComponentSwitch.hasLoosingBearing
    LoosingBearing = InitialParameter.LoosingBearing; hasLoosingBearing = true;
else
    hasLoosingBearing = false; 
end

if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasCustom') && InitialParameter.ComponentSwitch.hasCustom
    Custom = InitialParameter.Custom; hasCustom = true;
else
    hasCustom = false; 
end

% --- UPDATED: Speed-Dependent Bearing Flag ---
if isfield(InitialParameter, 'ComponentSwitch') && isfield(InitialParameter.ComponentSwitch, 'hasSpeedDependentBearing') && InitialParameter.ComponentSwitch.hasSpeedDependentBearing
    SpdBearing = InitialParameter.SpeedDependentBearing; hasSpdBearing = true;
else
    hasSpdBearing = false;
end

%% 1. Generate Key Points (Updated for Segmented Shafts)
keyPoints = cell(Shaft.amount,1);
for iShaft = 1:1:Shaft.amount
    % Add Segment Endpoints as Key Points
    segLengths = Shaft.segmentLength{iShaft};
    segmentBoundaries = [0; cumsum(segLengths(:))];
    
    % Initialize keyPoints with segment boundaries
    keyPoints{iShaft} = segmentBoundaries;
    
    % --- Add Component Locations ---
    % Record disk
    position = find(Disk.inShaftNo == iShaft);
    keyPointsDisk = Disk.positionOnShaftDistance(position);
    
    % Record ordinary bearing 
    position = find(Bearing.inShaftNo == iShaft);
    keyPointsBearing = Bearing.positionOnShaftDistance(position);
    
    % Record Speed-Dependent Bearing
    if hasSpdBearing
        position = find(SpdBearing.inShaftNo == iShaft);
        keyPointsSpdBearing = SpdBearing.positionOnShaftDistance(position);
    else
        keyPointsSpdBearing = [];
    end
    
    % Record rub
    if hasRub
        position = find(RubImpact.inShaftNo == iShaft);
        keyPointsRub = RubImpact.positionOnShaftDistance(position);
    else
        keyPointsRub = []; 
    end
    
    % Record intermediate bearing
    if hasInterBearing
        position = find(InterBearing.betweenShaftNo == iShaft);
        keyPointsInterBearing = InterBearing.positionOnShaftDistance(position);
    else 
        keyPointsInterBearing = []; 
    end
    
    % Record Customize function
    if hasCustom
        position = find(Custom.inShaftNo == iShaft);
        keyPointsCustom = Custom.positionOnShaftDistance(position);
    else 
        keyPointsCustom = []; 
    end
    
    % --- Assemble and Sort ---
    keyPoints{iShaft} = sort([  keyPoints{iShaft};...
                                keyPointsDisk;...
                                keyPointsRub;...
                                keyPointsInterBearing;...
                                keyPointsBearing; ...
                                keyPointsSpdBearing; ... 
                                keyPointsCustom]);
          
    % --- Merge Close Points ---
    ii = 1;
    totalL = Shaft.totalLength(iShaft);
    while ii < length(keyPoints{iShaft})
        distance = abs( keyPoints{iShaft}(ii) - keyPoints{iShaft}(ii+1) );
        isTooClose = distance < totalL/THRESHOLD_COEFFICIENT;
        if isTooClose
            keyPoints{iShaft}(ii+1) = [];
        else
            ii = ii + 1;
        end 
    end 
end 

%% 2. Check Meshing Parameters
if isAutoMesh
    switch gridFineness
        case 'low', FINENESS = 1;
        case 'middle', FINENESS = 4;
        case 'high', FINENESS = 10;
        otherwise, error('Please input: low, middle or high')
    end 
else
    for iShaft = 1:1:Shaft.amount
        isMatch = length(manualGrid{iShaft}) == ( length(keyPoints{iShaft}) - 1 );
        if ~isMatch
            error(['Manual grid dimension mismatch. Segments: ', num2str(length(keyPoints{iShaft}) - 1)]);
        end 
    end 
end 

%% 3. Perform Meshing (Generate Nodes)
rowSegmentNum = zeros(Shaft.amount,1);
nodeDistance = cell(Shaft.amount,1);
for iShaft = 1:1:Shaft.amount
    rowSegmentNum(iShaft) = length(keyPoints{iShaft}) - 1; 
    
    if isAutoMesh
        standardLength = Shaft.totalLength(iShaft)/ FINENESS;
        elementNum = ceil( diff(keyPoints{iShaft}) ./ standardLength );
    else
        elementNum = manualGrid{iShaft};
    end
    
    nodeDistance{iShaft} = [];
    for iSegment = 1:1:rowSegmentNum(iShaft)
        nodesInSegment = linspace( keyPoints{iShaft}(iSegment),...
                                   keyPoints{iShaft}(iSegment+1),...
                                   elementNum(iSegment)+1 );
        nodeDistance{iShaft} = [nodeDistance{iShaft}, nodesInSegment];
    end
    
    nodeDistance{iShaft} = unique(nodeDistance{iShaft});
end 

%% 4. Generate Node Structure
nodeNum = sum( cellfun(@length,nodeDistance) );
Node = struct('name', cell(nodeNum,1), 'onShaftNo', [], 'onShaftDistance', [],...
              'diskNo', [], 'bearingNo', [], 'isBearing', [], 'dof', []);

% Initialize optional fields
for iNode = 1:nodeNum
    if hasRub, Node(iNode).rubImpactNo = []; end
    if hasInterBearing, Node(iNode).interBearingNo = []; end
    if hasCoupling, Node(iNode).couplingNo = []; end
    if hasLoosingBearing, Node(iNode).isLoosingBearing = []; end
    if hasCustom, Node(iNode).CustomNo = []; end
    if hasSpdBearing, Node(iNode).speedDependentBearingNo = []; end 
end
 
iShaft = 1;
previousShaftNodeNum = 0;
for iNode = 1:1:nodeNum
    Node(iNode).name = iNode;
    Node(iNode).onShaftNo = iShaft;
    distanceHere = nodeDistance{iShaft}(iNode - previousShaftNodeNum);
    Node(iNode).onShaftDistance = distanceHere;
    Node(iNode).dof = Shaft.dofOfEachNodes(iShaft);
    
    [previousShaftNodeNum, iShaft] = judgeShaftEnd(previousShaftNodeNum,...
                                     nodeDistance, iShaft, iNode); 
end 

%% 5. Generate ShaftElement Structure
currentElemIdx = 0;
shaftElementNum = nodeNum-Shaft.amount;
ShaftElement(shaftElementNum) = struct();
for iShaft = 1:Shaft.amount
    nodes = nodeDistance{iShaft};
    numElemsInShaft = length(nodes) - 1;
    
    inSegLengths  = Shaft.segmentLength{iShaft};
    inOuterR      = Shaft.outerRadius{iShaft};
    inInnerR      = Shaft.innerRadius{iShaft};
    inOuterRStiff = Shaft.outerRadiusStiff{iShaft};
    inInnerRStiff = Shaft.innerRadiusStiff{iShaft};
    inDensity     = Shaft.density{iShaft};
    inE           = Shaft.elasticModulus{iShaft};
    inPoisson     = Shaft.poissonRatio{iShaft};
    inEcc         = Shaft.eccentricity{iShaft};
    inEccPhase    = Shaft.eccentricityPhase{iShaft};
    
    inSegBoundaries = [0; cumsum(inSegLengths(:))];
    
    nodeOffset = 0;
    for k = 1:(iShaft-1)
        nodeOffset = nodeOffset + length(nodeDistance{k});
    end
    
    for j = 1:numElemsInShaft
        currentElemIdx = currentElemIdx + 1;
        
        nodeL = nodeOffset + j;     
        nodeR = nodeOffset + j + 1; 
        
        zL = nodes(j);
        zR = nodes(j+1);
        L_elem = zR - zL;
        zMid = (zL + zR) / 2;
        
        segIdx = find(zMid >= inSegBoundaries(1:end-1) - 1e-9 & ...
                      zMid <= inSegBoundaries(2:end) + 1e-9, 1, 'first');
        
        if isempty(segIdx)
            error('Meshing Error: Element midpoint not found in any segment.');
        end
        
        ShaftElement(currentElemIdx).NodeNo = [nodeL, nodeR];
        ShaftElement(currentElemIdx).ShaftNo = iShaft;
        ShaftElement(currentElemIdx).Length = L_elem;
        
        ShaftElement(currentElemIdx).outerRadius = inOuterR(segIdx);
        ShaftElement(currentElemIdx).innerRadius = inInnerR(segIdx);
        ShaftElement(currentElemIdx).outerRadiusStiff = inOuterRStiff(segIdx);
        ShaftElement(currentElemIdx).innerRadiusStiff = inInnerRStiff(segIdx);
        ShaftElement(currentElemIdx).density = inDensity(segIdx);
        ShaftElement(currentElemIdx).elasticModulus = inE(segIdx);
        ShaftElement(currentElemIdx).poissonRatio = inPoisson(segIdx);
        ShaftElement(currentElemIdx).eccentricity = inEcc(segIdx);
        ShaftElement(currentElemIdx).eccentricityPhase = inEccPhase(segIdx);
    end
end

%% 6. Match Components to Nodes
[Disk,Node] = matchElement(Disk, 'disk', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft);
[Bearing,Node] = matchElement(Bearing, 'bearing', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft);
if hasRub, [RubImpact,Node] = matchElement(RubImpact, 'rubImpact', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft); end
if hasCoupling, [Coupling,Node] = matchElement(Coupling, 'couplingMisalignment', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft); end
if hasInterBearing, [InterBearing,Node] = matchElement(InterBearing, 'interBearing', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft); end
if hasCustom, [Custom,Node] = matchElement(Custom, 'custom', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft); end

% Match Speed-Dependent Bearing
if hasSpdBearing
    [SpdBearing,Node] = matchElement(SpdBearing, 'speedDependentBearing', Node, nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft); 
end

%% 7. Handle Loosing Bearing
if hasLoosingBearing
    loosingBearingIndex = LoosingBearing.inBearingNo;
    loosingNode = Bearing.positionOnShaftNode(loosingBearingIndex);
    [Node(loosingNode).isLoosingBearing] = deal(true);
    for iNode = 1:1:nodeNum
        if isempty(Node(iNode).isLoosingBearing)
            Node(iNode).isLoosingBearing = false;
        end
    end
end

%% 8. Add Bearing Nodes (Append to list)
[nodeNum, Bearing, Node] = addNode(nodeNum, Bearing, Node, hasLoosingBearing, 'bearing');

if hasInterBearing
    [nodeNum, InterBearing, Node] = addNode(nodeNum, InterBearing, Node, hasLoosingBearing, 'interBearing');
end

% Generate Internal Nodes for Speed-Dependent Bearing (if mass > 0)
if hasSpdBearing
    [nodeNum, SpdBearing, Node] = addNode(nodeNum, SpdBearing, Node, hasLoosingBearing, 'speedDependentBearing');
end

%% 9. DOF Calculation
dofOfEachNode = [Node.dof];
dofInterval = zeros(nodeNum,2);
for iNode = 1:1:nodeNum
    endDof = sum( dofOfEachNode(1:iNode) );
    startDof = endDof - dofOfEachNode(iNode) +1;
    dofInterval(iNode,:) = [startDof, endDof];
end

dofOnNodeNo = zeros(sum([Node.dof]),1);
dofNo = 1;
for iNode = 1:1:nodeNum
    for iDof = 1:1:Node(iNode).dof
        dofOnNodeNo(dofNo) = Node(iNode).name;
        dofNo = dofNo + 1;
    end
end

%% 10. Output Packaging
Parameter = InitialParameter;
% Mesh Output
Mesh.Node = Node;
Mesh.ShaftElement = ShaftElement; 
Mesh.keyPointsDistance = keyPoints;
Mesh.nodeDistance = nodeDistance;
Mesh.nodeNum = nodeNum;
Mesh.dofNum = sum([Node.dof]);
Mesh.dofOnNodeNo = dofOnNodeNo;
Mesh.dofInterval = dofInterval;
Parameter.Mesh = Mesh;

% Component Output
Parameter.Disk = Disk;
Parameter.Bearing = Bearing;
if hasRub, Parameter.RubImpact = RubImpact; end
if hasCoupling, Parameter.CouplingMisalignment = Coupling; end
if hasInterBearing, Parameter.IntermediateBearing = InterBearing; end
if hasCustom, Parameter.Custom = Custom; end

% Output Speed-Dependent Bearing
if hasSpdBearing, Parameter.SpeedDependentBearing = SpdBearing; end

end

%% Subfunctions
function [Element,Node] = matchElement(Element, elementName, Node,...
                          nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft)
    if strcmp(elementName,'interBearing'), columnNum = 2; else, columnNum = 1; end
    Element.positionOnShaftNode = zeros(Element.amount,columnNum);    
    shaftNo = 1; previousShaftNode = 0; iElement = 1;
    for nodeNo = 1:1:nodeNum
        distanceH = nodeDistance{shaftNo}(nodeNo - previousShaftNode);
        if strcmp(elementName,'interBearing'), isShaftHere = shaftNo == Element.betweenShaftNo;
        else, isShaftHere = shaftNo == Element.inShaftNo; end
        
        isDistanceHere = abs(distanceH - Element.positionOnShaftDistance) < Shaft.totalLength(shaftNo)/THRESHOLD_COEFFICIENT;
        isElementHere = isShaftHere & isDistanceHere; 
        isElementHereNum = length(isElementHere(isElementHere == true));
        
        if isElementHereNum == 1 
            indexElement = find(isElementHere == true);
            Element.positionOnShaftNode(indexElement) = nodeNo;
            switch elementName
                case 'disk', Node(nodeNo).diskNo = iElement;
                case 'bearing', Node(nodeNo).bearingNo = iElement;
                case 'rubImpact', Node(nodeNo).rubImpactNo = iElement;
                case 'couplingMisalignment', Node(nodeNo).couplingNo = iElement;
                case 'custom', Node(nodeNo).customNo = iElement;
                case 'speedDependentBearing', Node(nodeNo).speedDependentBearingNo = iElement; 
                case 'interBearing'
                    if sum(isElementHere(:,1)) == 1, Node(nodeNo).interBearingNo = iElement; 
                    else, iElement = iElement - 1; end 
            end 
            iElement = iElement +1;
        elseif isElementHereNum > 1 
            error('Distance between elements too close, adjust THRESHOLD_COEFFICIENT.');
        end 
        [previousShaftNode, shaftNo] = judgeShaftEnd(previousShaftNode, nodeDistance, shaftNo, nodeNo);
    end 
    if strcmp(elementName,'interBearing')
        index1Column = Element.positionOnShaftNode(:,1); index2Column = Element.positionOnShaftNode(:,2);
        for iInterBearing = 1:1:length(index1Column)
            Node(index2Column(iInterBearing)).interBearingNo = Node(index1Column(iInterBearing)).interBearingNo;
        end
    end 
end 

function [previousShaftNodeNum, iShaft] = judgeShaftEnd(previousShaftNodeNum, nodeDistance, iShaft, iNode)
    isShaftEnd = iNode == (previousShaftNodeNum + length(nodeDistance{iShaft}) );
    if isShaftEnd
        previousShaftNodeNum = previousShaftNodeNum + length(nodeDistance{iShaft});
        iShaft = iShaft +1;
    end 
end 

function [nodeNum, Element, Node] = addNode(nodeNum, Element, Node, hasLoosingBearing, elementType)
    indexHasMass = find(Element.mass' ~= 0); 
    [indexHasMassRow, indexHasMassCol] = find(Element.mass' ~= 0); 
    bearingNodeNum = length(indexHasMass);
    bearingMassColumnNum = size(Element.mass,2);
    Element.positionNode = zeros(Element.amount,bearingMassColumnNum); 
    for iBearingNode = 1:1:bearingNodeNum
        newNodeNo = nodeNum + iBearingNode;
        onNodeNo = Element.positionOnShaftNode(indexHasMassCol(iBearingNode),:);
        Node(newNodeNo).onShaftNo = [Node(onNodeNo).onShaftNo];
        Node(newNodeNo).onShaftDistance = [Node(onNodeNo).onShaftDistance];
        
        if strcmp(elementType, 'interBearing')
            Node(newNodeNo).interBearingNo = unique([Node(onNodeNo).interBearingNo]);
        elseif strcmp(elementType, 'speedDependentBearing')
            Node(newNodeNo).speedDependentBearingNo = Node(onNodeNo).speedDependentBearingNo;
        else
            Node(newNodeNo).bearingNo = Node(onNodeNo).bearingNo; 
        end
        
        if hasLoosingBearing, Node(newNodeNo).isLoosingBearing = Node(onNodeNo).isLoosingBearing; end
        Node(newNodeNo).name = newNodeNo; 
        Node(newNodeNo).dof = Element.dofOfEachNodes(indexHasMassCol(iBearingNode),indexHasMassRow(iBearingNode));
        Node(newNodeNo).isBearing = true;
        Element.positionNode(indexHasMassCol(iBearingNode),indexHasMassRow(iBearingNode)) = newNodeNo;
    end
    for iNodee = 1:1:nodeNum
        if isempty(Node(iNodee).isBearing), Node(iNodee).isBearing = false; end
    end
    nodeNum = nodeNum + bearingNodeNum;
end