%% femDisk - Generate FEM matrices for disk components in rotor systems
%
% This function assembles global mass, gyroscopic, and transient matrices 
% along with gravity vectors and unbalance data for disk elements in rotor 
% dynamics models.
%
%% Syntax
%  [M, G, N, Fg, DiskUnbalanceData] = femDisk(Disk, nodeDof)
%
%% Description
% |femDisk| constructs finite element matrices for disk components in 
% rotor systems. The function:
% * Computes disk mass and inertia properties
% * Generates mass and gyroscopic matrices
% * Creates gravity force vectors
% * Aggregates unbalance information (Magnitude/Phase) per node
% * Supports multiple disk configurations
%
%% Input Arguments
% * |Disk| - Disk properties structure:
%   * |amount|              % Number of disks (scalar)
%   * |dofOfEachNodes|      % DOF per mounted node [N×1 vector]
%   * |radius|              % Outer radii [m] [N×1 vector]
%   * |thickness|           % Axial thicknesses [m] [N×1 vector]
%   * |density|             % Material densities [kg/m³] [N×1 vector]
%   * |positionOnShaftNode| % Mounting node indices [N×1 vector]
%   * N: Number of disks
%
% * |nodeDof| - DOF counts per system node [M×1 vector], M = number of nodes
%
%% Output Arguments
% * |M|  % Global mass matrix [sparse n_total×n_total]
% * |G|  % Global gyroscopic matrix [sparse n_total×n_total]
% * |N|  % Nonlinear matrix [sparse n_total×n_total]
% * |Fg| % Gravity force vector [n_total×1]
% * |DiskUnbalanceData| % Unbalance data matrix [N_valid×3]:
%       Column 1: Node ID
%       Column 2: Magnitude [kg·m] (Mass × Eccentricity)
%       Column 3: Phase [rad]
%   * n_total: Total DOF of rotor system = sum(nodeDof)
%
%% Matrix Assembly Process
% 1. Element Generation:
%    * For each disk, calls |diskElement| to compute:
%      - Mass matrix (Me)
%      - Gyroscopic matrix (Ge)
%      - Transient matrix (Ne)
%      - Gravity force vector (Fge)
%      - Unbalance Information (Magnitude/Phase)
% 2. Global Matrix Initialization:
%    * Creates zero matrices of size sum(nodeDof)
% 3. Position Mapping:
%    * Determines DOF positions using |findIndex|
% 4. Assembly:
%    * Adds each disk's matrices to global positions via |addElementIn|
%    * Aggregates and filters unbalance data
%
%% Physical Modeling
% * Mass Calculation:
%   $ m = \pi \rho t (r_o^2 - r_i^2) $
% * Inertia Terms:
%   $ I_d = \frac{1}{12}m(3(r_o^2 + r_i^2) + t^2) $ (Diametral)
%   $ I_p = \frac{1}{2}m(r_o^2 + r_i^2) $ (Polar)
% * Matrix Structures:
%   * Mass matrix combines translational and rotational terms
%   * Gyroscopic matrix accounts for polar inertia effects
%
%% Implementation Notes
% * Eccentricity Handling:
%   * Returns unbalance data matrix [NodeID, Mag, Phase]
%   * Filters out negligible unbalance values (< 1e-15)
% * Position Mapping:
%   * Uses |findIndex| for DOF position calculation
% * Matrix Assembly:
%   * Utilizes |addElementIn| for efficient sparse matrix construction
%
%% Example
% % Configure disk parameters
% diskCfg = struct('amount', 2, ...
%                  'dofOfEachNodes', [4; 4], ...
%                  'outerRadius', [0.15; 0.12], ...
%                  'innerRadius', [0; 0], ...
%                  'thickness', [0.025; 0.02], ...
%                  'density', [7850; 7850], ...
%                  'eccentricity', [1e-3; 1e-3],...
%                  'positionOnShaftNode', [3; 5]);
% % System DOF configuration
% nodeDOF = [4,4,4,4,4,4,4,4,4,4]'; 
% % Generate disk matrices
% [M, G, ~, Fg, UnbalData] = femDisk(diskCfg, nodeDOF);
%
%% See Also
% diskElement, addElementIn, findIndex, femShaft, femBearing
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%



function [M, G, N, Fg, DiskUnbalanceData] = femDisk(Disk, nodeDof)


% Initialize element matrix storage
Me = cell(Disk.amount,1); 
Ge = cell(Disk.amount,1);
Ne = cell(Disk.amount,1);
Fge = cell(Disk.amount,1);

% Initialize Unbalance Data Storage
% Format: [NodeID, Magnitude, Phase]
rawUnbalance = zeros(Disk.amount, 3);

% Prepare temporary struct for getStructPiece (Legacy support)
Temporary = rmfield(Disk, 'amount'); 

%% 1. Element Generation Loop
for iDisk = 1:1:Disk.amount
    % Extract single disk data
    ADisk = getStructPiece(Temporary, iDisk, []);
    
    % Generate matrices and unbalance info for this disk
    [Me{iDisk}, Ge{iDisk}, Ne{iDisk}, Fge{iDisk}, unbalInfo] = diskElement(ADisk); 
    
    % Store Unbalance Data
    % Disk.positionOnShaftNode contains the Global Node ID
    nodeID = Disk.positionOnShaftNode(iDisk);
    rawUnbalance(iDisk, :) = [nodeID, unbalInfo.Magnitude, unbalInfo.Phase];
end

%% 2. Global Matrix Assembly
dofNum = sum(nodeDof);
M = zeros(dofNum, dofNum);
G = zeros(dofNum, dofNum);
N = zeros(dofNum, dofNum);
Fg = zeros(dofNum, 1);

% Calculate the position of disk element in global matrix
% (Maps Node ID to Matrix Indices)
diskOnDofPosition = findIndex(Disk.positionOnShaftNode, nodeDof);

for iDisk = 1:1:Disk.amount
   M = addElementIn(M, Me{iDisk}, diskOnDofPosition(iDisk, :));
   G = addElementIn(G, Ge{iDisk}, diskOnDofPosition(iDisk, :));
   N = addElementIn(N, Ne{iDisk}, diskOnDofPosition(iDisk, :));
   Fg = addElementIn(Fg, Fge{iDisk}, [diskOnDofPosition(iDisk,1), 1]);
end

%% 3. Process Disk Unbalance Data
% Filter out disks with zero unbalance (Optional, but good for consistency)
tolerance = 1e-15;
validIndices = rawUnbalance(:, 2) > tolerance;

% Output formatted matrix
DiskUnbalanceData = rawUnbalance(validIndices, :);

end