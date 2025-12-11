%% bearingElementInter - Generate stiffness/damping matrices for non-mass intermediate bearings
%
% This function constructs partitioned stiffness and damping matrices for 
% intermediate bearings without concentrated mass connecting two shaft nodes.
% It supports full 2x2 matrices allowing for cross-coupling effects.
%
%% Syntax
%  [Ke, Ce] = bearingElementInter(ABearing)
%
%% Description
% |bearingElementInter| calculates partitioned stiffness and damping matrices 
% for massless intermediate bearings connecting two shaft nodes. The function:
% * Models anisotropic bearing properties including cross-coupling (k_xy, etc.)
% * Creates partitioned matrices for node-node connectivity (FEM assembly)
% * Handles different DOF counts on connected shaft nodes
%
%% Input Arguments
% * |ABearing| - Bearing properties structure:
%   * |dofOnShaftNode|    % DOF counts on connected shaft nodes [1×2 vector]
%   * |stiffness|         % Horizontal stiffness (kV) [N/m]
%   * |stiffnessVertical| % Vertical stiffness (kW) [N/m]
%   * |stiffnessHV|       % (Optional) Cross-stiffness Force_H/Disp_V [N/m]
%   * |stiffnessVH|       % (Optional) Cross-stiffness Force_V/Disp_H [N/m]
%   * |damping|           % Horizontal damping (cV) [N·s/m]
%   * |dampingVertical|   % Vertical damping (cW) [N·s/m]
%   * |dampingHV|         % (Optional) Cross-damping Force_H/Vel_V [N·s/m]
%   * |dampingVH|         % (Optional) Cross-damping Force_V/Vel_H [N·s/m]
%
%% Output Arguments
% * |Ke| - Partitioned stiffness matrix (2×2 cell array):
%   * |Ke{1,1}| - Stiffness matrix for start node DOF (Self)
%   * |Ke{1,2}| - Coupling stiffness matrix (Start -> End)
%   * |Ke{2,1}| - Coupling stiffness matrix (End -> Start)
%   * |Ke{2,2}| - Stiffness matrix for end node DOF (Self)
% * |Ce| - Partitioned damping matrix (2×2 cell array with same structure)
%
%% Matrix Construction Rules
% 1. Local Element Matrix (Kin/Cin):
%   * Full 2x2 matrix: |Kin = [kV, kHV; kVH, kW]|
% 2. Partitioning (Relative Displacement Formulation):
%   * |K11| =  |Kin| (Node 1 resists motion)
%   * |K12| = |-Kin| (Node 1 pulled by Node 2)
%   * |K21| = |-Kin| (Node 2 pulled by Node 1)
%   * |K22| =  |Kin| (Node 2 resists motion)
%
%% Example
% % Configure intermediate bearing with cross-coupling
% bearingCfg.dofOnShaftNode = [4, 4];
% bearingCfg.stiffness = 1e6;            % k_xx
% bearingCfg.stiffnessVertical = 1.2e6;  % k_yy
% bearingCfg.stiffnessHV = 2e5;          % k_xy (New)
% bearingCfg.stiffnessVH = -2e5;         % k_yx (New)
% bearingCfg.damping = 500;
% bearingCfg.dampingVertical = 600;
% [Ke, Ce] = bearingElementInter(bearingCfg);
%
%% See Also
% bearingElement, bearingElementInterMass, addElementIn
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function [Ke, Ce] = bearingElementInter(ABearing)
%% Check input stiffness and damping
% Filters out zero values for main diagonal terms validation
ABearing.stiffness(find(ABearing.stiffness==0)) = [];
ABearing.stiffnessVertical(find(ABearing.stiffnessVertical==0)) = [];
ABearing.damping(find(ABearing.damping==0)) = [];
ABearing.dampingVertical(find(ABearing.dampingVertical==0)) = [];

% Validate Stiffness
if isempty(ABearing.stiffness)
    ABearing.stiffness = 0;
else
    if length(ABearing.stiffness) ~= 1
        error('too much input stiffness for bearing without mass')
    end
end

if isempty(ABearing.stiffnessVertical)
    ABearing.stiffnessVertical = 0;
else
    if length(ABearing.stiffnessVertical) ~= 1
        error('too much input vertical stiffness for bearing without mass')
    end
end

% Validate Damping
if isempty(ABearing.damping)
    ABearing.damping = 0;
else
    if length(ABearing.damping) ~= 1
        error('too much input damping for bearing without mass')
    end
end

if isempty(ABearing.dampingVertical)
    ABearing.dampingVertical = 0;
else
    if length(ABearing.dampingVertical) ~= 1
        error('too much input vertical damping for bearing without mass')
    end
end

%% Constants Extraction
kV = ABearing.stiffness;
cV = ABearing.damping;
kW = ABearing.stiffnessVertical;
cW = ABearing.dampingVertical;

% Cross-coupling terms (Check existence, default to 0)
if isfield(ABearing, 'stiffnessHV'), kHV = ABearing.stiffnessHV(1); else, kHV = 0; end
if isfield(ABearing, 'stiffnessVH'), kVH = ABearing.stiffnessVH(1); else, kVH = 0; end
if isfield(ABearing, 'dampingHV'),   cHV = ABearing.dampingHV(1);   else, cHV = 0; end
if isfield(ABearing, 'dampingVH'),   cVH = ABearing.dampingVH(1);   else, cVH = 0; end

dof1 = ABearing.dofOnShaftNode(1);
dof2 = ABearing.dofOnShaftNode(2);

%% Stiffness Matrix Construction
% Construct the full 2x2 local stiffness matrix including cross-terms
Kin = [ kV,  kHV; ...
        kVH, kW ];
 
% Initialize blocks
K11 = zeros(dof1);      K12 = zeros(dof1, dof2);
K21 = K12';             K22 = zeros(dof2);

% Assembly based on relative displacement: F = K * (x1 - x2)
% Node 1 forces
K11 = addElementIn(K11,  Kin, [1,1]);  
K12 = addElementIn(K12, -Kin, [1,1]);
% Node 2 forces
K21 = addElementIn(K21, -Kin, [1,1]); 
K22 = addElementIn(K22,  Kin, [1,1]);

Ke = {K11, K12;...
      K21, K22 };
  
%% Damping Matrix Construction
% Construct the full 2x2 local damping matrix including cross-terms
Cin = [ cV,  cHV; ...
        cVH, cW ];
 
% Initialize blocks
C11 = zeros(dof1);      C12 = zeros(dof1, dof2);
C21 = C12';             C22 = zeros(dof2);

% Assembly based on relative velocity
C11 = addElementIn(C11,  Cin, [1,1]);  
C12 = addElementIn(C12, -Cin, [1,1]);
C21 = addElementIn(C21, -Cin, [1,1]); 
C22 = addElementIn(C22,  Cin, [1,1]);

Ce = {C11, C12;...
      C21, C22 };
end