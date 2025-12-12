%% standardizeInputParameter - Standardize input parameters for backward compatibility
%
% This function ensures that the input parameter structure conforms to the 
% latest version (V2) specification, enabling the software to process both 
% legacy (V1) and modern input formats seamlessly.
%
%% Syntax
%  NewParam = standardizeInputParameter(OldParam)
%
%% Description
% |standardizeInputParameter| checks the format of the input structure. 
% If the input is in the legacy format (V1, based on double arrays for 
% shaft properties), it converts it to the modern format (V2, based on 
% cell arrays for segmented shafts).
%
% It automatically handles:
% * Conversion of vector data to cell array data.
% * Initialization of new fields: |outerRadiusStiff|, |innerRadiusStiff|, 
%   |eccentricity|, and |eccentricityPhase|.
% * Setting default values for new fields (e.g., Stiffness Radius defaults 
%   to Geometric Radius).
%
%% Input Arguments
% * |OldParam| - Structure containing system parameters (Status, Shaft, Disk, etc.)
%   Can be either Version 1 (Legacy) or Version 2 (Segment-based).
%
%% Output Arguments
% * |NewParam| - Standardized structure with V2 formatting.
%   * Shaft.outerRadius      - Cell array format
%   * Shaft.outerRadiusStiff - Populated (new field)
%   * Shaft.eccentricity     - Populated (new field)
%
%% See Also
%  inputEssentialParameterTwinSpool, inputEssentialParameterTwinSpool2
%
% Copyright (c) 2021-2025 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%
%%
function NewParam = standardizeInputParameter(OldParam)

NewParam = OldParam;
Shaft = NewParam.Shaft;

% Check if conversion is needed
% V1: outerRadius is a double vector (e.g., [0.01; 0.02])
% V2: outerRadius is a cell array (e.g., {[...]; [...]})
if iscell(Shaft.outerRadius)
    % Already V2 format, do nothing
    return;
end

% fprintf('Detected Legacy Input Format. Converting to Segment-based Format...\n');

%% Perform Conversion from V1 to V2
numShafts = Shaft.amount;

% Initialize new Cell Arrays
new_segmentLength     = cell(numShafts, 1);
new_outerRadius       = cell(numShafts, 1);
new_innerRadius       = cell(numShafts, 1);
new_outerRadiusStiff  = cell(numShafts, 1);
new_innerRadiusStiff  = cell(numShafts, 1);
new_density           = cell(numShafts, 1);
new_elasticModulus    = cell(numShafts, 1);
new_poissonRatio      = cell(numShafts, 1);
new_eccentricity      = cell(numShafts, 1);
new_eccentricityPhase = cell(numShafts, 1);

for i = 1:numShafts
    % --- Geometry & Material ---
    % In V1, the whole shaft is 1 segment.
    % So segmentLength = totalLength
    new_segmentLength{i}  = Shaft.totalLength(i);
    
    % V1 had scalar/vector properties per shaft, V2 needs a vector inside a cell
    new_outerRadius{i}    = Shaft.outerRadius(i);
    new_innerRadius{i}    = Shaft.innerRadius(i);
    new_density{i}        = Shaft.density(i);
    new_elasticModulus{i} = Shaft.elasticModulus(i);
    new_poissonRatio{i}   = Shaft.poissonRatio(i);
    
    % --- New Fields in V2 (Handling Backward Compatibility) ---
    
    % 1. Stiffness Radii: Default to Geometric Radii
    % Assumption: If not specified (old code), Stiffness = Physical Geometry
    new_outerRadiusStiff{i} = Shaft.outerRadius(i);
    new_innerRadiusStiff{i} = Shaft.innerRadius(i);
    
    % 2. Shaft Distributed Eccentricity: Default to 0
    % Assumption: Old models treat shafts as perfectly balanced lines
    new_eccentricity{i}      = 0;
    new_eccentricityPhase{i} = 0;
end

% Assign back to structure
NewParam.Shaft.segmentLength     = new_segmentLength;
NewParam.Shaft.outerRadius       = new_outerRadius;
NewParam.Shaft.innerRadius       = new_innerRadius;
NewParam.Shaft.outerRadiusStiff  = new_outerRadiusStiff; % New field
NewParam.Shaft.innerRadiusStiff  = new_innerRadiusStiff; % New field
NewParam.Shaft.density           = new_density;
NewParam.Shaft.elasticModulus    = new_elasticModulus;
NewParam.Shaft.poissonRatio      = new_poissonRatio;
NewParam.Shaft.eccentricity      = new_eccentricity;         % New field
NewParam.Shaft.eccentricityPhase = new_eccentricityPhase; % New field

% Note: totalLength, dofOfEachNodes, rayleighDamping are compatible 
% between versions (still simple vectors/scalars), so we leave them unchanged.

end