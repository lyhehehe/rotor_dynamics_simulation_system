%% bearingElement - Generate stiffness/damping matrices for massless bearings
%
% This function constructs local stiffness and damping matrices for bearing 
% elements without concentrated mass in rotor dynamics models. It supports 
% full 2x2 stiffness/damping matrices allowing for cross-coupling effects.
%
%% Syntax
%  [Ke, Ce] = bearingElement(ANBearing)
%
%% Description
% |bearingElement| calculates stiffness and damping matrices for massless 
% bearing elements. The function:
% * Validates bearing property inputs (Specific zero-filtering logic)
% * Constructs anisotropic stiffness/damping matrices with cross-coupling
% * Automatically expands matrices to match nodal DOF count
%
%% Input Arguments
% * |ANBearing| - Bearing properties structure with fields:
%   * |dofOnShaftNode|    % DOF count per bearing node (scalar integer)
%   * |stiffness|         % Horizontal stiffness (kV) [N/m]
%   * |stiffnessVertical| % Vertical stiffness (kW) [N/m]
%   * |stiffnessHV|       % (Optional) Cross-stiffness (Force_H / Disp_V) [N/m]
%   * |stiffnessVH|       % (Optional) Cross-stiffness (Force_V / Disp_H) [N/m]
%   * |damping|           % Horizontal damping (cV) [N·s/m]
%   * |dampingVertical|   % Vertical damping (cW) [N·s/m]
%   * |dampingHV|         % (Optional) Cross-damping (Force_H / Vel_V) [N·s/m]
%   * |dampingVH|         % (Optional) Cross-damping (Force_V / Vel_H) [N·s/m]
%
%   Note: If cross-coupling fields are omitted, they default to zero.
%
%% Output Arguments
% * |Ke| - Local stiffness matrix (n×n), where n = dofOnShaftNode
% * |Ce| - Local damping matrix (n×n), where n = dofOnShaftNode
%
%% Matrix Construction Rules
% 1. Anisotropic modeling:
%    * Supports full 2x2 element matrices: [k_V, k_HV; k_VH, k_W]
% 2. Matrix expansion:
%    * Stiffness/damping terms placed in the first 2 DOF positions (H, V)
%    * Remaining DOF filled with zeros
% 3. Zero value handling:
%    * Specific filtering applied to main stiffness/damping inputs
%
%% Example
% % Create massless bearing parameters
% bearing = struct('dofOnShaftNode', 4, ...
%                  'stiffness', 1e8, ...
%                  'stiffnessVertical', 1.2e8, ...
%                  'stiffnessHV', 2e6, ...   % Cross-term
%                  'damping', 500, ...
%                  'dampingVertical', 600);
% % Generate stiffness and damping matrices
% [Ke, Ce] = bearingElement(bearing);
%
%% See Also
% bearingElementMass, diskElement, shaftElement, assemblyGlobalMatrix
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
function [Ke, Ce] = bearingElement(ANBearing)
%% check input stiffness and damping (Original Logic Preserved)
ANBearing.stiffness(find(ANBearing.stiffness==0)) = [];
ANBearing.stiffnessVertical(find(ANBearing.stiffnessVertical==0)) = [];
ANBearing.damping(find(ANBearing.damping==0)) = [];
ANBearing.dampingVertical(find(ANBearing.dampingVertical==0)) = [];

if isempty(ANBearing.stiffness)
    ANBearing.stiffness = 0;
else
    % check length
    if length(ANBearing.stiffness) ~= 1
        error('too much input stiffness for bearing without mass')
    end
end
if isempty(ANBearing.damping)
    ANBearing.damping = 0;
else
    % check length
    if length(ANBearing.damping) ~= 1
        error('too much input damping for bearing without mass')
    end
end

%% constants
% Main diagonal terms
kV = ANBearing.stiffness;         % V direction: horizontal
kW = ANBearing.stiffnessVertical; % W direction: vertical
cV = ANBearing.damping;
cW = ANBearing.dampingVertical;

% Cross-coupling terms (New added features)
% Check if fields exist, otherwise default to zeros
if isfield(ANBearing, 'stiffnessHV'), kHV = ANBearing.stiffnessHV(1); else, kHV = 0; end
if isfield(ANBearing, 'stiffnessVH'), kVH = ANBearing.stiffnessVH(1); else, kVH = 0; end
if isfield(ANBearing, 'dampingHV'),   cHV = ANBearing.dampingHV(1);   else, cHV = 0; end
if isfield(ANBearing, 'dampingVH'),   cVH = ANBearing.dampingVH(1);   else, cVH = 0; end

dof = ANBearing.dofOnShaftNode;

%% generate stiffness matrix of bearing element
% Modified to include cross-coupling terms
Ke = [ kV,  kHV;
       kVH, kW ];
   
% Expand stiffness matrix (fill remaining DOFs with zeros)
% Note: Using dof - 2 because the bearing matrix is explicitly 2x2
Ke = blkdiag( Ke, zeros(dof) ); 

%% generate damping matrix of bearing element
% Modified to include cross-coupling terms
Ce = [ cV,  cHV;
       cVH, cW ];

% Expand damping matrix (fill remaining DOFs with zeros)
Ce = blkdiag( Ce, zeros(dof) ); 

end