%% bearingElementInterMass - Generate FEM matrices for mass-bearing elements with multi-DOF
%
% This function constructs partitioned FEM matrices for intermediate bearings 
% with multiple concentrated masses and cross-coupling stiffness/damping.
%
%% Syntax
%  [Me, Ke, Ce, Fge] = bearingElementInterMass(AMBearing)
%
%% Description
% |bearingElementInterMass| calculates partitioned mass, stiffness, damping, 
% and gravity force matrices for intermediate bearings with:
% * Multiple concentrated masses between two shaft segments (Inner/Outer)
% * Support for anisotropic properties and cross-coupling (k_xy, c_xy, etc.)
%
%% Input Arguments
% * |AMBearing| - Bearing configuration structure:
%   * |dofOnShaftNode|    % DOF counts for connected shaft nodes [1×2 double] (Inner, Outer)
%   * |dofOfEachNodes|    % DOF counts per bearing mass node [1×n double]
%   * |mass|              % Concentrated masses [kg] 
%   * |stiffness|         % Horizontal stiffness (kV) [1×(n+1)]
%   * |stiffnessVertical| % Vertical stiffness (kW) [1×(n+1)]
%   * |stiffnessHV|       % (Optional) Cross-stiffness Force_H/Disp_V [1×(n+1)]
%   * |stiffnessVH|       % (Optional) Cross-stiffness Force_V/Disp_H [1×(n+1)]
%   * |damping|           % Horizontal damping (cV)
%   * |dampingVertical|   % Vertical damping (cW)
%   * |dampingHV|         % (Optional) Cross-damping
%   * |dampingVH|         % (Optional) Cross-damping
%
%% Output Arguments
% * |Ke| - Partitioned stiffness matrix (1×7 cell array):
%   * |Ke{1}| - Inner shaft stiffness
%   * |Ke{2}| - Inner shaft to first mass coupling
%   * |Ke{3}| - First mass to inner shaft coupling
%   * |Ke{4}| - Outer shaft stiffness
%   * |Ke{5}| - Outer shaft to last mass coupling
%   * |Ke{6}| - Last mass to outer shaft coupling
%   * |Ke{7}| - Mass chain stiffness (inter-mass connections)
% * |Ce| - Partitioned damping matrix (1×7 cell array with same structure as Ke)
% * |Fge| - Gravity force vector for mass nodes [n×1 double]
%
%% Matrix Construction Rules
% 1. Mass handling:
%   * Non-zero masses are filtered from input
%   * Masses are placed on diagonal positions in |Me{2,2}|
% 2. Stiffness/damping chains:
%   * Components connect sequential mass nodes
%   * Horizontal/vertical terms remain uncoupled
%   * Terminal connections linked to shaft nodes
% 3. Gravity forces:
%   * Applied to vertical DOF of each mass node
%   * Magnitude: -9.8 * mass (downward direction)
%
%% Dimension Requirements
% * |mass|, |dofOfEachNodes|, |stiffness|, |damping| must have compatible lengths:
%   * |n = length(nonzero_mass)|
%   * |stiffness| and |damping| must have length n+1
%   * |dofOfEachNodes| must have length ≥ n
%
%% Example
% % Configure intermediate bearing with two masses
% bearingConfig = struct('dofOnShaftNode', [4,4], ...
%                       'mass', [3.5, 2.1], ...      % Two non-zero masses
%                       'stiffness', [1e6, 8e5, 7e5], ... % Three stiffness components
%                       'damping', [1e3, 8e2, 7e2], ...
%                       'dampingVertical', [1e3, 8e2, 7e2], ...
%                       'stiffnessVertical', [1e6, 8e5, 7e5], ...
%                       'dofOfEachNodes', [2,2]);   % Two mass nodes with 2 DOF each
% % Generate partitioned matrices
% [Me, Ke, Ce, Fge] = bearingElementInterMass(bearingConfig);
%
%% See Also
% bearingElementInter, addElementIn
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.

function [Me, Ke, Ce, Fge] = bearingElementInterMass(AMBearing)
%% 1. Initial & Data Extraction
m = AMBearing.mass;
% Main diagonal
kV = AMBearing.stiffness;
kW = AMBearing.stiffnessVertical;
cV = AMBearing.damping;
cW = AMBearing.dampingVertical;
% Cross-coupling (Check existence, default to 0)
if isfield(AMBearing, 'stiffnessHV'), kHV = AMBearing.stiffnessHV; else, kHV = zeros(size(kV)); end
if isfield(AMBearing, 'stiffnessVH'), kVH = AMBearing.stiffnessVH; else, kVH = zeros(size(kV)); end
if isfield(AMBearing, 'dampingHV'),   cHV = AMBearing.dampingHV;   else, cHV = zeros(size(cV)); end
if isfield(AMBearing, 'dampingVH'),   cVH = AMBearing.dampingVH;   else, cVH = zeros(size(cV)); end

dofShaft = AMBearing.dofOnShaftNode;   % [dof_S1, dof_S2]
dofBearing = AMBearing.dofOfEachNodes; 

m = m(m~=0); 
massNum = length(m);

% Ensure property vector lengths are (massNum + 1)
expectedLen = massNum + 1;
if length(kV) ~= expectedLen
    kV = kV(1:expectedLen); kW = kW(1:expectedLen);
    cV = cV(1:expectedLen); cW = cW(1:expectedLen);
    kHV = kHV(1:expectedLen); kVH = kVH(1:expectedLen);
    cHV = cHV(1:expectedLen); cVH = cVH(1:expectedLen);
    dofBearing = dofBearing(1:massNum);
end

% Total DOF definition
% Structure of Global Matrix K: [Shaft1_DOFs, Shaft2_DOFs, Mass1_DOFs, ..., MassN_DOFs]
dofShaftTotal = sum(dofShaft);
dofBearingTotal = sum(dofBearing);
dofTotal = dofShaftTotal + dofBearingTotal;

K = zeros(dofTotal, dofTotal);
C = zeros(dofTotal, dofTotal);
Fge = zeros(dofBearingTotal, 1);

%% 2. Global Matrix Assembly
% Concept: We iterate through ELEMENTS (Springs/Dampers).
% Total Elements = massNum + 1
% Element 1: Shaft1 <-> Mass1
% Element i: Mass(i-1) <-> Mass(i)
% Element End: Mass(n) <-> Shaft2

% Define indices for Shafts in the global matrix
idxShaft1 = 1:dofShaft(1);
idxShaft2 = (dofShaft(1) + 1) : dofShaftTotal;
idxMassStart = dofShaftTotal; % Offset for mass nodes

% Pre-calculate mass node indices
massIndices = cell(1, massNum);
currentIdx = idxMassStart;
for i = 1:massNum
    massIndices{i} = currentIdx + (1:dofBearing(i));
    currentIdx = currentIdx + dofBearing(i);
end

for iElem = 1:(massNum + 1)
    
    % Get 2x2 element matrices (Stiffness & Damping)
    [Ke_elem, Ce_elem] = getElemMat(iElem);
    
    % --- Identify Node indices (Left and Right of the element) ---
    if iElem == 1
        % 1st Element: Connects Shaft 1 (Left) -> Mass 1 (Right)
        idxL = idxShaft1;        
        idxR = massIndices{1};   
    elseif iElem == (massNum + 1)
        % Last Element: Connects Mass N (Left) -> Shaft 2 (Right)
        idxL = massIndices{end}; 
        idxR = idxShaft2;        
    else
        % Intermediate Element: Connects Mass i-1 (Left) -> Mass i (Right)
        idxL = massIndices{iElem - 1};
        idxR = massIndices{iElem};
    end
    
    % --- Assembly (Standard FEM) ---
    % 1. Add to diagonal (Self-stiffness)
    % Left Node
    K = addElementIn(K, Ke_elem, [idxL(1), idxL(1)]);
    C = addElementIn(C, Ce_elem, [idxL(1), idxL(1)]);
    % Right Node
    K = addElementIn(K, Ke_elem, [idxR(1), idxR(1)]);
    C = addElementIn(C, Ce_elem, [idxR(1), idxR(1)]);
    
    % 2. Add to off-diagonal (Coupling) -> Negative Sign
    % K_LR (Upper Right)
    K = addElementIn(K, -Ke_elem, [idxL(1), idxR(1)]);
    C = addElementIn(C, -Ce_elem, [idxL(1), idxR(1)]);
    % K_RL (Lower Left)
    K = addElementIn(K, -Ke_elem, [idxR(1), idxL(1)]);
    C = addElementIn(C, -Ce_elem, [idxR(1), idxL(1)]);
    
end

%% 3. Partitioning (Output Generation)
% The logic below relies on K being organized as [S1, S2, M1, M2...].
% We use mat2cell to slice it based on DOF counts.

% Prepare cell division vector: [dofS1, dofS2, dofM1, dofM2, ...]
splitVector = [dofShaft, dofBearing]; 

KeTemp = mat2cell(K, splitVector, splitVector);
CeTemp = mat2cell(C, splitVector, splitVector);

% Extract specific blocks based on function specification
% Indices in KeTemp: 1=S1, 2=S2, 3=M1, 4=M2 ... (Note: 2 is Outer Shaft)

if massNum == 1
    % Structure: 1=S1, 2=S2, 3=M1
    % Ke{1}=K_S1_S1, Ke{2}=K_S1_M1, Ke{3}=K_M1_S1
    % Ke{4}=K_S2_S2, Ke{5}=K_S2_M1, Ke{6}=K_M1_S2, Ke{7}=K_M1_M1
    Ke = {KeTemp{1,1}, KeTemp{1,3}, KeTemp{3,1}, ...
          KeTemp{2,2}, KeTemp{2,3}, KeTemp{3,2}, KeTemp{3,3}};
    Ce = {CeTemp{1,1}, CeTemp{1,3}, CeTemp{3,1}, ...
          CeTemp{2,2}, CeTemp{2,3}, CeTemp{3,2}, CeTemp{3,3}};
      
elseif massNum >= 2
    % Structure: 1=S1, 2=S2, 3=M1, ..., Last=Mn
    idxM_First = 3;
    idxM_Last  = 2 + massNum;
    
    % Mass Chain Block (Ke{7}): Need to extract the whole mass-mass block
    % This is K(dofShaftTotal+1:end, dofShaftTotal+1:end)
    K_MassChain = K(dofShaftTotal+1:end, dofShaftTotal+1:end);
    C_MassChain = C(dofShaftTotal+1:end, dofShaftTotal+1:end);
    
    Ke = {KeTemp{1,1}, KeTemp{1,idxM_First}, KeTemp{idxM_First,1}, ... % S1 related
          KeTemp{2,2}, KeTemp{2,idxM_Last},  KeTemp{idxM_Last,2},  ... % S2 related
          K_MassChain};                                                % Mass Chain
    Ce = {CeTemp{1,1}, CeTemp{1,idxM_First}, CeTemp{idxM_First,1}, ...
          CeTemp{2,2}, CeTemp{2,idxM_Last},  CeTemp{idxM_Last,2},  ...
          C_MassChain};
end

%% 4. Mass Matrix
Min = zeros(dofBearingTotal);
currentIdx = 0;
for im = 1:massNum
    Mi = eye(2) * m(im); % Diagonal mass
    Min = addElementIn(Min, Mi, [currentIdx+1, currentIdx+1]);
    currentIdx = currentIdx + dofBearing(im);
end
% Partitioned Output
M11 = zeros(dofShaftTotal);      
M12 = zeros(dofShaftTotal, dofBearingTotal);
M21 = M12';                    
M22 = Min; % Mass matrix for bearing nodes
Me = {M11, M12; M21, M22};

%% 5. Gravity
Fge = zeros(dofBearingTotal, 1);
for im = 1:massNum
    FgeHere = -9.8 * m(im);
    idxBase = sum(dofBearing(1:im-1));
    Fge(idxBase + 2) = FgeHere; % Apply to Vertical DOF (Local 2)
end

    %% Sub-function: Get Element Matrix (Includes Cross-Coupling)
    function [Ke_local, Ce_local] = getElemMat(idx)
        Ke_local = [kV(idx),  kHV(idx); ...
                    kVH(idx), kW(idx)];
        Ce_local = [cV(idx),  cHV(idx); ...
                    cVH(idx), cW(idx)];
    end

end