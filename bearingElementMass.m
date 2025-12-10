%% bearingElementMass - Generate FEM matrices for bearings with concentrated masses and cross-coupling
%
% This function constructs partitioned mass, stiffness, and damping matrices 
% along with gravity force vectors for bearing elements with concentrated masses 
% in rotor dynamics systems. It supports full 2x2 stiffness/damping matrices 
% allowing for cross-coupling effects (e.g., oil-film bearings).
%
%% Syntax
%  [Me, Ke, Ce, Fge] = bearingElementMass(AMBearing)
%
%% Description
% |bearingElementMass| calculates local matrices for bearing elements with 
% concentrated masses. The function:
% * Models anisotropic bearing properties including cross-coupling terms
%   (stiffness/damping coefficients k_xy, k_yx, etc.)
% * Handles multi-node mass configurations
% * Generates partitioned matrices for efficient global assembly
% * Computes gravity forces for mass components
%
%% Input Arguments
% * |AMBearing| - Bearing properties structure:
%   * |dofOfEachNodes|    % DOF counts per bearing node [1×K vector]
%   * |stiffness|         % Horizontal stiffness (kV: Force_H / Disp_H) [N/m] [1×M]
%   * |stiffnessVertical| % Vertical stiffness (kW: Force_V / Disp_V) [N/m] [1×M]
%   * |stiffnessHV|       % (Optional) Cross-stiffness (Force_H / Disp_V) [N/m] [1×M]
%   * |stiffnessVH|       % (Optional) Cross-stiffness (Force_V / Disp_H) [N/m] [1×M]
%   * |damping|           % Horizontal damping (cV) [N·s/m] [1×M]
%   * |dampingVertical|   % Vertical damping (cW) [N·s/m] [1×M]
%   * |dampingHV|         % (Optional) Cross-damping (Force_H / Vel_V) [N·s/m] [1×M]
%   * |dampingVH|         % (Optional) Cross-damping (Force_V / Vel_H) [N·s/m] [1×M]
%   * |mass|              % Concentrated masses [kg] [1×K vector of non-zero masses]
%   * |dofOnShaftNode|    % DOF count at shaft connection node [scalar]
%
%   Note: M = K + 1 (properties count = mass nodes + 1). 
%   If cross-coupling fields are omitted, they default to zero.
%
%% Output Arguments
% * |Me| - Partitioned mass matrix (2×2 cell array):
%   * |Me{1,1}|, |Me{1,2}|, |Me{2,1}| - Zero matrices (Shaft coupling)
%   * |Me{2,2}| - Lumped mass matrix for bearing nodes (Diagonal)
% * |Ke| - Partitioned stiffness matrix (2×2 cell array):
%   * |Ke{1,1}| - Shaft-shaft stiffness coupling (includes 1st bearing element)
%   * |Ke{1,2}| - Shaft-mass stiffness coupling
%   * |Ke{2,1}| - Mass-shaft stiffness coupling
%   * |Ke{2,2}| - Mass-mass stiffness coupling (includes cross-coupling terms)
% * |Ce| - Partitioned damping matrix (2×2 cell array with same structure as Ke)
% * |Fge| - Gravity force vector [N] (vector for bearing DOF, V-direction only)
%
%% Matrix Construction Rules
% 1. Mass handling:
%   * Non-zero masses are placed on diagonal positions in |Me{2,2}|
% 2. Stiffness/damping:
%   * Supports full 2x2 element matrices: [k_V, k_HV; k_VH, k_W]
%   * Elements are assembled in series: Shaft <-> Node1 <-> ... <-> NodeN <-> Ground
%   * Terminal connections (Element 1) linked to shaft node
%   * Final connections (Element N+1) linked to ground
% 3. Gravity forces:
%   * Applied only to vertical DOF of each mass node (-9.8 * mass)
%
%% Example
% % Configure bearing with cross-coupling (e.g., simplified oil film)
% bearingProps = struct('dofOfEachNodes', [2, 2], ...
%                      'stiffness',         [1e8, 5e7, 3e7], ... % k_xx
%                      'stiffnessVertical', [1.2e8, 6e7, 3.5e7], ... % k_yy
%                      'stiffnessHV',       [2e6, 1e6, 0.5e6], ...   % k_xy (New)
%                      'stiffnessVH',       [-2e6, -1e6, -0.5e6], ...% k_yx (New)
%                      'damping',           [500, 300, 200], ...
%                      'dampingVertical',   [600, 350, 250], ...
%                      'mass',              [3.5, 2.1], ...
%                      'dofOnShaftNode',    4);
% % Generate partitioned matrices
% [Me, Ke, Ce, Fg] = bearingElementMass(bearingProps);
%
%% See Also
% bearingElementInterMass, addElementIn, assembleLinear
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%


function [Me, Ke, Ce, Fge] = bearingElementMass(AMBearing)
%% 1. Initial & Data Extraction
% Constants extraction
m = AMBearing.mass;
% Main diagonal terms (V=Horizontal, W=Vertical)
kV = AMBearing.stiffness; 
kW = AMBearing.stiffnessVertical;
cV = AMBearing.damping;
cW = AMBearing.dampingVertical;

% Cross-coupling terms (New added features)
% Check if fields exist, otherwise default to zeros for backward compatibility
if isfield(AMBearing, 'stiffnessHV'), kHV = AMBearing.stiffnessHV; else, kHV = zeros(size(kV)); end
if isfield(AMBearing, 'stiffnessVH'), kVH = AMBearing.stiffnessVH; else, kVH = zeros(size(kV)); end
if isfield(AMBearing, 'dampingHV'),   cHV = AMBearing.dampingHV;   else, cHV = zeros(size(cV)); end
if isfield(AMBearing, 'dampingVH'),   cVH = AMBearing.dampingVH;   else, cVH = zeros(size(cV)); end

dofShaft = AMBearing.dofOnShaftNode;
dofBearingNodes = AMBearing.dofOfEachNodes;

% Filter non-zero masses
nonZeroIndices = m ~= 0;
m = m(nonZeroIndices); 
massNum = length(m);

% If input arrays are longer than needed, truncate them to match (massNum + 1 elements)
if length(AMBearing.stiffness) ~= (massNum + 1)
    % We need N+1 elements for N masses
    numElem = massNum + 1;
    kV = kV(1:numElem);   kW = kW(1:numElem);
    cV = cV(1:numElem);   cW = cW(1:numElem);
    kHV = kHV(1:numElem); kVH = kVH(1:numElem);
    cHV = cHV(1:numElem); cVH = cVH(1:numElem);
    dofBearingNodes = dofBearingNodes(1:massNum);
end

dofBearingTotal = sum(dofBearingNodes);
dofTotal = dofShaft + dofBearingTotal;

% Initialize Global Matrices
K = zeros(dofTotal, dofTotal);
C = zeros(dofTotal, dofTotal);
Fge = zeros(dofBearingTotal, 1);

%% 2. Matrix Assembly Loop
% Iterate through each mass node to assemble the system
% Connectivity: Shaft <-> Mass 1 <-> Mass 2 ... <-> Mass N <-> Ground

accumulatedDof = dofShaft; % Track the starting DOF index for the current bearing node

for im = 1:massNum
    
    % Current Node Info
    currentMassDof = dofBearingNodes(im);
    idxNode = accumulatedDof + (1:currentMassDof); % Indices of current node in global matrix
    
    % --- Element Identification ---
    % Element L: Connects "Previous Node" (or Shaft) to "Current Node"
    idxElemL = im; 
    [Ke_L, Ce_L] = getElemMat(idxElemL);
    
    % Element R: Connects "Current Node" to "Next Node" (or Ground)
    idxElemR = im + 1;
    [Ke_R, Ce_R] = getElemMat(idxElemR);
    
    % --- Assembly Logic ---
    
    % A. Handle Connection to Left (Previous Node or Shaft)
    if im == 1
        % Case: First Mass Node connecting to Shaft
        
        % 1. Shaft-Shaft Coupling (Top-Left of Global K)
        % Adds contribution of Element 1 to the shaft
        K = addElementIn(K, Ke_L, [1, 1]); 
        C = addElementIn(C, Ce_L, [1, 1]);
        
        % 2. Shaft-Node Coupling (Off-diagonal)
        % Note: Left coupling is -K_elem
        K = addElementIn(K, -Ke_L, [1, accumulatedDof+1]); % K_12
        K = addElementIn(K, -Ke_L, [accumulatedDof+1, 1]); % K_21
        C = addElementIn(C, -Ce_L, [1, accumulatedDof+1]); 
        C = addElementIn(C, -Ce_L, [accumulatedDof+1, 1]);
        
    else
        % Case: Intermediate Mass Node
        % Previous Node is Mass(im-1)
        prevDofCount = dofBearingNodes(im-1);
        idxPrev = (accumulatedDof - prevDofCount) + (1:prevDofCount);
        
        % Coupling with previous mass (Standard FEM: -K)
        % Note: The positive K for the previous node was added in the previous loop iteration
        K = addElementIn(K, -Ke_L, [accumulatedDof+1, idxPrev(1)]); % Lower-Left block
        C = addElementIn(C, -Ce_L, [accumulatedDof+1, idxPrev(1)]);
        
        % Symmetric coupling (Upper-Right block) handled by addElementIn if implemented symmetrically,
        % or we explicitly add the upper triangle part:
        K = addElementIn(K, -Ke_L, [idxPrev(1), accumulatedDof+1]); 
        C = addElementIn(C, -Ce_L, [idxPrev(1), accumulatedDof+1]);
    end
    
    % B. Handle Diagonal of Current Node (Self-Stiffness)
    % Contribution from Left Element + Right Element
    % K_node = K_L + K_R
    K = addElementIn(K, Ke_L + Ke_R, [idxNode(1), idxNode(1)]);
    C = addElementIn(C, Ce_L + Ce_R, [idxNode(1), idxNode(1)]);
    
    % C. Gravity Force
    % Gravity only acts on the Vertical direction (local index 2)
    % F = m * g (downward)
    FgeHere = -9.8 * m(im);
    % Map to the specific DOF in the bearing-only vector
    % Calculate local index in Fge vector
    idxLocalBase = sum(dofBearingNodes(1:im-1));
    Fge(idxLocalBase + 2) = FgeHere; 
    
    % Update DOF counter for next iteration
    accumulatedDof = accumulatedDof + currentMassDof;
    
end

%% 3. Partition Matrices (Output Formatting)
% Me, Ke, Ce structure:
% { Shaft-Shaft,  Shaft-Bearing }
% { Bearing-Shaft, Bearing-Bearing }

Ke = mat2cell(K, [dofShaft, dofBearingTotal], [dofShaft, dofBearingTotal]);
Ce = mat2cell(C, [dofShaft, dofBearingTotal], [dofShaft, dofBearingTotal]);

%% 4. Mass Matrix Assembly
% Mass is purely diagonal and lumped at nodes
M_bearing = zeros(dofBearingTotal, dofBearingTotal);
current_idx = 0;
for im = 1:massNum
    % 2x2 Mass block: diag([m, m])
    Mi = eye(2) * m(im);
    M_bearing = addElementIn(M_bearing, Mi, [current_idx+1, current_idx+1]);
    current_idx = current_idx + dofBearingNodes(im);
end

% Partitioned Mass Matrix (Shaft has no mass contribution in this element function)
M11 = zeros(dofShaft, dofShaft);
M12 = zeros(dofShaft, dofBearingTotal);
M21 = zeros(dofBearingTotal, dofShaft);
M22 = M_bearing;

Me = {M11, M12; M21, M22};


    %% Sub-function: Get Element Matrix
    % Retrieves the 2x2 Stiffness and Damping matrices for the i-th element
    % Includes full cross-coupling terms
    function [Ke_local, Ce_local] = getElemMat(idx)
        % Stiffness
        % [ kV   kHV ]
        % [ kVH  kW  ]
        Ke_local = [kV(idx),  kHV(idx); ...
                    kVH(idx), kW(idx)];
                    
        % Damping
        % [ cV   cHV ]
        % [ cVH  cW  ]
        Ce_local = [cV(idx),  cHV(idx); ...
                    cVH(idx), cW(idx)];
    end

end