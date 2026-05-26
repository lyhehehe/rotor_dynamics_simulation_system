%% updateSpeedDependentMatrices - Build global stiffness and damping matrices for speed-dependent bearings
%
% This function interpolates per-bearing stiffness and damping lookup tables
% at the current shaft speeds and assembles the contributions into the global
% sparse matrices for use in the frequency-domain solver.
%
%% Syntax
%  [K_out, C_out] = updateSpeedDependentMatrices(SpdBearing, K_base, C_base, currentSpeeds, dofInterval)
%
%% Description
% |updateSpeedDependentMatrices| loops over each speed-dependent bearing,
% determines its shaft operating speed from |currentSpeeds|, performs 1-D
% linear interpolation (clamped at the lookup table boundaries), and adds
% the resulting 2×2 [Kxx Kxy; Kyx Kyy] and [Cxx Cxy; Cyx Cyy] blocks to
% the global matrices. The function:
% * Returns the base matrices unchanged if |SpdBearing| is empty or has zero bearings
% * Issues a warning if any bearing has non-zero mass (mass matrix is not updated)
% * Uses manual interpolation — does not require the Curve Fitting Toolbox
%
%% Input Arguments
% * |SpdBearing|    - SpeedDependentBearing structure (from |inputSpeedDependentBearing*|)
% * |K_base|        - Base global stiffness matrix [sparse, dofNum×dofNum]
% * |C_base|        - Base global damping matrix [sparse, dofNum×dofNum]
% * |currentSpeeds| - Current operating speeds for each shaft [rad/s] [shaftNum×1 double]
% * |dofInterval|   - DOF index range for each node [N×2 integer array]
%
%% Output Arguments
% * |K_out| - Updated global stiffness matrix [sparse, dofNum×dofNum]
% * |C_out| - Updated global damping matrix [sparse, dofNum×dofNum]
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function [K_out, C_out] = updateSpeedDependentMatrices(SpdBearing, K_base, C_base, currentSpeeds, dofInterval)

K_out = K_base;
C_out = C_base;

if isempty(SpdBearing) || SpdBearing.amount == 0
    return;
end

if any(sum(SpdBearing.mass, 2) ~= 0)
    warning('SpeedDependentBearing contains non-zero mass. The global M matrix is NOT updated.');
end

speedVec = SpdBearing.speed;
numSpeeds = length(speedVec);

% Loop through each bearing
for iBearing = 1:SpdBearing.amount
    
    % --- CORRECTED: Determine physical operating speed for THIS specific bearing ---
    shaftIdx = SpdBearing.inShaftNo(iBearing);
    opSpeed = currentSpeeds(shaftIdx); % Operating speed of the shaft it's attached to
    
    % Manual 1D Linear Interpolation Logic with Clamp Boundary
    if opSpeed <= speedVec(1)
        idx1 = 1; idx2 = 1; alpha = 0;
    elseif opSpeed >= speedVec(end)
        idx1 = numSpeeds; idx2 = numSpeeds; alpha = 0;
    else
        idx1 = find(speedVec <= opSpeed, 1, 'last');
        idx2 = idx1 + 1;
        alpha = (opSpeed - speedVec(idx1)) / (speedVec(idx2) - speedVec(idx1));
    end
    
    % --- Interpolation ---
    % Stiffness
    kxx = SpdBearing.stiffness(iBearing, idx1) + alpha * (SpdBearing.stiffness(iBearing, idx2) - SpdBearing.stiffness(iBearing, idx1));
    kyy = SpdBearing.stiffnessVertical(iBearing, idx1) + alpha * (SpdBearing.stiffnessVertical(iBearing, idx2) - SpdBearing.stiffnessVertical(iBearing, idx1));
    kxy = SpdBearing.stiffnessHV(iBearing, idx1) + alpha * (SpdBearing.stiffnessHV(iBearing, idx2) - SpdBearing.stiffnessHV(iBearing, idx1));
    kyx = SpdBearing.stiffnessVH(iBearing, idx1) + alpha * (SpdBearing.stiffnessVH(iBearing, idx2) - SpdBearing.stiffnessVH(iBearing, idx1));
    
    % Damping
    cxx = SpdBearing.damping(iBearing, idx1) + alpha * (SpdBearing.damping(iBearing, idx2) - SpdBearing.damping(iBearing, idx1));
    cyy = SpdBearing.dampingVertical(iBearing, idx1) + alpha * (SpdBearing.dampingVertical(iBearing, idx2) - SpdBearing.dampingVertical(iBearing, idx1));
    cxy = SpdBearing.dampingHV(iBearing, idx1) + alpha * (SpdBearing.dampingHV(iBearing, idx2) - SpdBearing.dampingHV(iBearing, idx1));
    cyx = SpdBearing.dampingVH(iBearing, idx1) + alpha * (SpdBearing.dampingVH(iBearing, idx2) - SpdBearing.dampingVH(iBearing, idx1));
    
    % Matrices
    Ke = [kxx, kxy; kyx, kyy];
    Ce = [cxx, cxy; cyx, cyy];
          
    % Direct assembly
    nodeNo = SpdBearing.positionOnShaftNode(iBearing);
    startDof = dofInterval(nodeNo, 1);
    dofIdx = [startDof, startDof + 1];
    
    K_out(dofIdx, dofIdx) = K_out(dofIdx, dofIdx) + Ke;
    C_out(dofIdx, dofIdx) = C_out(dofIdx, dofIdx) + Ce;
end

end