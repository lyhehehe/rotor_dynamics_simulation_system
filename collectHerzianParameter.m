%% collectHerzianParameter - Collect essential parameters for Hertzian contact force calculation
%
% This function extracts and organizes geometric and mechanical parameters 
% from both standard bearings and intermediate bearings to support nonlinear 
% Hertzian contact force computations.
%
%% Syntax
%  HerzianParameter = collectHerzianParameter(Mesh, ComponentSwitch, Shaft, Bearing)
%  HerzianParameter = collectHerzianParameter(Mesh, ComponentSwitch, Shaft, Bearing, InterBearing)
%
%% Description
% |collectHerzianParameter| processes system configuration data to generate 
% a unified parameter structure for contact mechanics. The function:
% * Identifies bearings requiring Hertzian force calculation
% * Maps Degrees of Freedom (DOFs) for inner and outer races
% * Handles ground-to-shaft and shaft-to-shaft (intermediate) connections
% * Consolidates geometric properties (radii, clearance) and material constants
% * Automates shaft index mapping for relative rotational speed calculations
%
%% Input Arguments
% * |Mesh| - Discretization data structure containing:
%   * |dofNum|: Total number of DOFs
%   * |dofInterval|: Node-to-DOF mapping matrix
% * |ComponentSwitch| - Logic flags for system components:
%   * |hasIntermediateBearing|: Flag for shaft-to-shaft bearings
% * |Shaft| - Shaft configuration containing |amount| (total shaft count)
% * |Bearing| - Standard bearing data structure:
%   * |isHertzian|: Logical array for Hertzian model activation
%   * |radiusInnerRace|, |radiusOuterRace|: Race radii [m]
%   * |rollerNum|: Number of rolling elements
%   * |clearance|: Radial clearance [m]
%   * |contactStiffness|: Hertzian contact stiffness [N/m^n]
%   * |coefficient|: Load-deflection exponent (e.g., 1.5 for ball bearings)
% * |InterBearing| - Intermediate bearing data structure (optional):
%   * Follows similar fields as |Bearing| with additional shaft-linking data
%
%% Output Arguments
% * |HerzianParameter| - Consolidated structure for force calculation:
%   * |nb|: Number of rollers [1 × N]
%   * |ri|, |ro|: Inner/outer race radii [1 × N]
%   * |delta0|: Radial clearances [1 × N]
%   * |kHertz|: Contact stiffness constants [1 × N]
%   * |n|: Contact exponents [1 × N]
%   * |omegaiNo|, |omegaoNo|: Shaft indices for speed referencing [1 × N]
%   * |hertzDof|: DOF indices for inner/outer race displacement [N × 2]
%   * |hertzianNum|: Total count of active Hertzian contact pairs
%
%% Parameter Mapping Logic
% 1. **Normal Bearings**: 
%    * Maps inner race to the specified shaft node.
%    * Maps outer race to ground (|dofNum + 1|) or a bearing housing mass node.
% 2. **Intermediate Bearings**: 
%    * Maps both races to respective nodes on the two connected shafts.
%    * Determines relative rotation by tracking |innerShaftNo| vs |outerShaftNo|.
%
%% See Also
% calculateResponse, dynamicEquation, hertzianForce
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
function HerzianParameter = collectHerzianParameter(Mesh, ComponentSwitch, Shaft, Bearing, InterBearing)

arguments
    Mesh 
    ComponentSwitch 
    Shaft 
    Bearing 
    InterBearing = []
end

% load some parameters
dofNum = Mesh.dofNum;
dofInterval = Mesh.dofInterval;
shaftNum = Shaft.amount;


%% calculate essential parameters for normal bearing
% the number of the hertzian force on the normal bearing
hertzianNumN = sum(Bearing.isHertzian); 

% initial
omegaiNoN = zeros(hertzianNumN, 1); % saving the inner shaft No for Herzian force
omegaoNoN = zeros(hertzianNumN, 1); % saving the outer shaft No for Herzian force

% hertzDof: is an index n*2 matrix. The first column saves the dof No. of
% the Inner shaft, the second column saves the dof No. of the Outer shaft.
% If this functioin is used to calculte the Hertzian force of a Bearing
% connecting the ground, the second column would expect a "dofNum+1". The
% i-th bearing info should be saved in the i-th row.
hertzDofN = zeros(hertzianNumN, 2); % hertzDof for normal bearing

riN = zeros(hertzianNumN, 1);
roN = zeros(hertzianNumN, 1);
nbN = zeros(hertzianNumN, 1);
delta0N = zeros(hertzianNumN, 1);
kHertzN = zeros(hertzianNumN, 1);
nN = zeros(hertzianNumN, 1);

% calculate
iHertz = 1;
for iBearing = 1:1:Bearing.amount
    if Bearing.isHertzian(iBearing)
        % for inner shaft
        shaftNo = Bearing.inShaftNo(iBearing);
        omegaiNoN(iHertz) = shaftNo;
        nodeNo = Bearing.positionOnShaftNode(iBearing);
        hertzDofN(iHertz, 1) = dofInterval(nodeNo, 1); % the x-direction dof of the node connecting this bearing
        
        % for outer shaft
        omegaoNoN(iHertz) = shaftNum + 1; % indicates the rotation speed of outer race equals to 0 (see hertzianForce.m)
        hasMass = sum(Bearing.mass(iBearing, :)) ~= 0;
        if hasMass
            massNodeNo = Bearing.positionNode(iBearing, 1);
            hertzDofN(iHertz, 2) = dofInterval(massNodeNo, 1);
        else
            hertzDofN(iHertz, 2) = dofNum + 1; % indicates the outer race connecting the ground (see hertzianForce.m)
        end % end if hasMass
        
        % for usual parameters
        riN(iHertz) = Bearing.radiusInnerRace(iBearing);
        roN(iHertz) = Bearing.radiusOuterRace(iBearing);
        nbN(iHertz) = Bearing.rollerNum(iBearing);
        delta0N(iHertz) = Bearing.clearance(iBearing);
        kHertzN(iHertz) = Bearing.contactStiffness(iBearing);
        nN(iHertz) = Bearing.coefficient(iBearing);
        
        % update
        iHertz = iHertz + 1;
    end % end if
end % end for


%% calcualte essential parameters for intermediate bearing
if ComponentSwitch.hasIntermediateBearing
    % the number of the hertzian force on the intermediate bearing
    hertzianNumI = sum(InterBearing.isHertzian); 
    
    % initial
    omegaiNoI = zeros(hertzianNumI, 1);
    omegaoNoI = zeros(hertzianNumI, 1);
    hertzDofI = zeros(hertzianNumI, 2); % hertzDof for intermediate bearing
    riI = zeros(hertzianNumI, 1);
    roI = zeros(hertzianNumI, 1);
    nbI = zeros(hertzianNumI, 1);
    delta0I = zeros(hertzianNumI, 1);
    kHertzI = zeros(hertzianNumI, 1);
    nI = zeros(hertzianNumI, 1);
    
    % calculate
    % generate outer shaft No based on the InterBearing.innershaftNo
    outerShaftNo = zeros(InterBearing.amount, 1); % initial
    for iBearing = 1:1:InterBearing.amount
        if InterBearing.innerShaftNo(iBearing)~=InterBearing.betweenShaftNo(iBearing, 1)
            outerShaftNo(iBearing) = InterBearing.betweenShaftNo(iBearing, 1);
        else
            outerShaftNo(iBearing) = InterBearing.betweenShaftNo(iBearing, 2);
        end % end if
    end % end for
    
    % calculate omegaiNoI, omegaoNoI, hertzDofI
    iHertz = 1;
    for iBearing = 1:1:InterBearing.amount
        if InterBearing.isHertzian(iBearing)
            % for inner shaft
            shaftNo = InterBearing.innerShaftNo(iBearing); % inner shaft No
            omegaiNoI(iHertz) = shaftNo;
            nodeNoIndex = find(InterBearing.betweenShaftNo(iBearing,:)==shaftNo);
            nodeNo = InterBearing.positionOnShaftNode(iBearing, nodeNoIndex);
            hertzDofI(iHertz, 1) = dofInterval(nodeNo, 1); % the x-direction dof of the node connecting this bearing
            
            % for outer shaft
            omegaoNoI(iHertz) = outerShaftNo(iBearing);
            hasMass = sum(InterBearing.mass(iBearing, :)) ~= 0;
            if hasMass
                if InterBearing.isHertzianTop(iBearing)
                    massNodeNo = InterBearing.positionNode(iBearing, 1);
                else
                    massNodeNoIndex = find(InterBearing.positionNode(iBearing, :),1,'last');
                    massNodeNo = InterBearing.positionNode(iBearing, massNodeNoIndex);
                end % end if top
                hertzDofI(iHertz, 2) = dofInterval(massNodeNo, 1);
            else
                nodeNoIndex = find(InterBearing.betweenShaftNo(iBearing,:)==outerShaftNo(iBearing));
                nodeNo = InterBearing.positionOnShaftNode(iBearing, nodeNoIndex);
                hertzDofI(iHertz, 2) = dofInterval(nodeNo, 1);
            end % end if hasMass
            
            % for usual parameters
            riI(iHertz) = InterBearing.radiusInnerRace(iBearing);
            roI(iHertz) = InterBearing.radiusOuterRace(iBearing);
            nbI(iHertz) = InterBearing.rollerNum(iBearing);
            delta0I(iHertz) = InterBearing.clearance(iBearing);
            kHertzI(iHertz) = InterBearing.contactStiffness(iBearing);
            nI(iHertz) = InterBearing.coefficient(iBearing);
            
            % update
            iHertz = iHertz + 1; 
        end % end if
    end % end for
end % end if (intermediate bearing)


%% Combining parameters in both Bearing and Intermediate Bearing 
if ComponentSwitch.hasIntermediateBearing
    nb = [nbN; nbI]';
    ri = [riN; riI]';
    ro = [roN; roI]';
    delta0 = [delta0N; delta0I]';
    kHertz = [kHertzN; kHertzI]';
    n = [nN; nI]';
    omegaiNo = [omegaiNoN; omegaiNoI]';
    omegaoNo = [omegaoNoN; omegaoNoI]';
    hertzDof = [hertzDofN; hertzDofI]';
    hertzianNum = hertzianNumN + hertzianNumI;
else
    nb = nbN';
    ri = riN';
    ro = roN';
    delta0 = delta0N';
    kHertz = kHertzN';
    n = nN';
    omegaiNo = omegaiNoN';
    omegaoNo = omegaoNoN';
    hertzDof = hertzDofN';
    hertzianNum = hertzianNumN;
end

% output Herzian parameters
HerzianParameter.nb = nb;
HerzianParameter.ri = ri;
HerzianParameter.ro = ro;
HerzianParameter.delta0 = delta0;
HerzianParameter.kHertz = kHertz;
HerzianParameter.n = n;
HerzianParameter.omegaiNo = omegaiNo;
HerzianParameter.omegaoNo = omegaoNo;
HerzianParameter.hertzDof = hertzDof';
HerzianParameter.hertzianNum = hertzianNum;

end