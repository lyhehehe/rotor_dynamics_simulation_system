%% calculateCampbell - Calculate Campbell diagram and critical speeds for a multi-shaft rotor system
%
% This function computes natural frequencies across a speed sweep and identifies
% critical speeds by intersecting mode frequency curves with engine-order
% excitation lines.
%
%% Syntax
%  [newEig, criticalSpeeds] = calculateCampbell(Parameter, exciteRad)
%  [newEig, criticalSpeeds] = calculateCampbell(Parameter, exciteRad, Name=Value)
%
%% Description
% |calculateCampbell| solves the generalized eigenvalue problem at each
% speed in |exciteRad| to track natural frequencies across the operating range.
% The function:
% * Builds the first-order state-space form and extracts eigenvalues
% * Optionally applies mode tracking using MAC or slope-based filtering
% * Plots the Campbell diagram and marks critical speeds if requested
% * Supports optional inclusion of the gyroscopic matrix
%
%% Input Arguments
% * |Parameter| - System configuration structure containing:
%   * |Matrix|: Global matrices (mass, stiffness, damping, gyroscopic)
%   * |Mesh|: Mesh data including |dofNum| and |dofInterval|
%   * |Status|: Speed ratio |ratio| for multi-shaft gyroscopic scaling
%   * |Shaft|: Shaft properties including |amount|
% * |exciteRad| - Speed sweep vector [rad/s] [1×N double]
%
%% Name-Value Arguments
% * |isPlot|          - Plot Campbell diagram (default: false) [logical]
% * |isFilter|        - Enable mode tracking across speed steps (default: false) [logical]
% * |filterMethod|    - Mode tracking algorithm: 'MAC' or 'slope' (default: 'MAC') [char]
% * |tolRad|          - Rigid-body mode rejection threshold [rad/s] (default: 0.1) [double]
% * |isUseGyroMatrix| - Include gyroscopic matrix in eigenvalue problem (default: true) [logical]
%
%% Output Arguments
% * |newEig|         - Natural frequency matrix [rad/s] [nModes×N double]
% * |criticalSpeeds| - Sorted critical speed vector [rad/s] where mode lines cross excitation lines [1×K double]
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function [newEig, criticalSpeeds] = calculateCampbell(Parameter, exciteRad, options)

    arguments
        Parameter struct
        exciteRad (1,:) double
        options.isPlot (1,1) logical = false
        options.isFilter (1,1) logical = false
        options.filterMethod (1,:) char {mustBeMember(options.filterMethod, {'slope', 'MAC'})} = 'MAC'
        options.tolRad (1,1) double = 1e-1
        options.isUseGyroMatrix (1,1) logical = true  % added option
    end
    
    % 1. Extract matrices and initialization
    M = Parameter.Matrix.mass;
    C = Parameter.Matrix.damping;
    K = Parameter.Matrix.stiffness;
    G_base = Parameter.Matrix.gyroscopic;
    
    dofNum = Parameter.Mesh.dofNum;
    eigNum = 2 * dofNum;
    exciteRadNum = length(exciteRad);
    speedRatio = Parameter.Status.ratio;
    shaftNum = Parameter.Shaft.amount;
    
    % Get shaft DOF range
    Node = Parameter.Mesh.Node;
    dofInterval = Parameter.Mesh.dofInterval;
    shaftDof = zeros(shaftNum, 2);
    for iShaft = 1:shaftNum
        IShaftNode = Node([Node.onShaftNo] == iShaft & [Node.isBearing] == false);
        startNode = min([IShaftNode.name]);
        endNode = max([IShaftNode.name]);
        shaftDof(iShaft,:) = [dofInterval(startNode,1), dofInterval(endNode,2)];
    end
    
    % Pre-allocate arrays
    rawEigVals = zeros(eigNum, exciteRadNum);
    if options.isUseGyroMatrix && options.isFilter && strcmp(options.filterMethod, 'MAC')
        rawEigVecs = cell(1, exciteRadNum);
    end
    
    % Precompute matrix divisions outside the loop
    invM_K = M \ K;
    invM_C = M \ C;
    
    % 2. Main Calculation 
    if options.isUseGyroMatrix
        for iRad = 1:exciteRadNum
            iBasicSpeed = exciteRad(iRad);
            G = G_base; 
            
            % Update G with speed ratio for each shaft
            for iShaft = 1:shaftNum
                if iShaft == 1
                    shaftSpeed = iBasicSpeed;
                else
                    shaftSpeed = iBasicSpeed * speedRatio(iShaft-1);
                end
                rng = shaftDof(iShaft,1) : shaftDof(iShaft,2);
                G(rng, rng) = shaftSpeed * G(rng, rng);
            end
            
            % Calculate matrix division for G
            invM_G = M \ G;
            
            % Assemble State-Space Matrix A (Using original sign convention: C-G)
            A = [-(invM_C - invM_G), -invM_K; ...
                 eye(dofNum),        zeros(dofNum, dofNum)];
            
            % Calculate eigenvalues (and eigenvectors if MAC is needed)
            if options.isFilter && strcmp(options.filterMethod, 'MAC')
                [V, D] = eig(full(A));
                rawEigVals(:, iRad) = diag(D);
                rawEigVecs{iRad} = V;
            else
                rawEigVals(:, iRad) = eig(full(A));
            end
        end
        
    else
        % No Gyroscopic effects: A is constant, compute only once
        A = [-invM_C,     -invM_K; ...
             eye(dofNum), zeros(dofNum, dofNum)];
             
        singleRawVals = eig(full(A));
        rawEigVals = repmat(singleRawVals, 1, exciteRadNum); % Duplicate across all speeds
    end
    
    % Extract natural frequencies (absolute value of imaginary parts, rad/s)
    freqs = abs(imag(rawEigVals));
    
    % 3. Mode Sorting / Tracking (Filtering)
    trackedEig = zeros(size(freqs));
    
    if ~options.isUseGyroMatrix
        % No tracking needed for horizontal lines, just sort numerically
        trackedEig = sort(freqs, 1);
        
    else
        % Original tracking logic for gyro-enabled calculation
        if ~options.isFilter
            trackedEig = sort(freqs, 1);
            
        elseif strcmp(options.filterMethod, 'MAC')
            trackedEig(:, 1) = freqs(:, 1);
            [trackedEig(:, 1), sortIdx] = sort(trackedEig(:, 1));
            prev_V = rawEigVecs{1}(:, sortIdx);
            
            for iRad = 2:exciteRadNum
                curr_freqs = freqs(:, iRad);
                curr_V = rawEigVecs{iRad};
                
                MAC = zeros(eigNum, eigNum);
                for j = 1:eigNum
                    for k = 1:eigNum
                        MAC(j,k) = abs(prev_V(:,j)' * curr_V(:,k))^2 / ...
                                   (norm(prev_V(:,j))^2 * norm(curr_V(:,k))^2);
                    end
                end
                
                matchIdx = zeros(eigNum, 1);
                for j = 1:eigNum
                    [~, maxIdx] = max(MAC(j, :));
                    matchIdx(j) = maxIdx;
                    MAC(:, maxIdx) = -1; 
                end
                
                trackedEig(:, iRad) = curr_freqs(matchIdx);
                prev_V = curr_V(:, matchIdx);
            end
            
        elseif strcmp(options.filterMethod, 'slope')
            trackedEig(:, 1) = sort(freqs(:, 1));
            
            if exciteRadNum > 1
                trackedEig(:, 2) = sort(freqs(:, 2));
                step = exciteRad(2) - exciteRad(1);
                prevSlope = (trackedEig(:, 2) - trackedEig(:, 1)) / step;
                
                for iRad = 3:exciteRadNum
                    step = exciteRad(iRad) - exciteRad(iRad-1);
                    predictedFreqs = trackedEig(:, iRad-1) + prevSlope * step;
                    curr_freqs = freqs(:, iRad);
                    
                    matchIdx = zeros(eigNum, 1);
                    used = false(eigNum, 1);
                    
                    for j = 1:eigNum
                        diffs = abs(curr_freqs - predictedFreqs(j));
                        diffs(used) = inf; 
                        [~, minIdx] = min(diffs);
                        matchIdx(j) = minIdx;
                        used(minIdx) = true;
                    end
                    
                    trackedEig(:, iRad) = curr_freqs(matchIdx);
                    prevSlope = (trackedEig(:, iRad) - trackedEig(:, iRad-1)) / step;
                end
            end
        end
    end
    
    % 4. Final Output Formatting & Fake Mode Deletion
    % Since eigenvalues are complex conjugates, every pair is duplicated in magnitude.
    % We take every 2nd row to return N unique modes instead of 2N.
    newEig = trackedEig(1:2:end, :);
    tolRad = options.tolRad; 
    
    if ~options.isUseGyroMatrix
        % Simplified filtering for constant horizontal lines
        validModeIdx = newEig(:, 1) > tolRad;
        newEig = newEig(validModeIdx, :);
        
    else
        % Filter out numerical rigid body modes with zero stiffness using y-intercept extrapolation
        if exciteRadNum > 1
            dOmega = exciteRad(2) - exciteRad(1);
            dFreq = newEig(:, 2) - newEig(:, 1);
            initialSlope = dFreq / dOmega;
            
            % Extrapolate the initial slope to Omega = 0 to find the static natural frequency intercept
            staticFreqIntercept = newEig(:, 1) - initialSlope * exciteRad(1);
            
            validModeIdx = abs(staticFreqIntercept) > tolRad;
            newEig = newEig(validModeIdx, :);
            
        elseif exciteRadNum == 1 && exciteRad(1) == 0
            validModeIdx = newEig(:, 1) > tolRad;
            newEig = newEig(validModeIdx, :);
        end
    end
    
    % 4.5. Calculate critical speeds (intersections between mode lines and synchronous excitation lines)
    criticalPoints = []; % Internal array to store [x, y] coordinate pairs
    speedRatio_full = [1; speedRatio(:)];
    
    for iShaft = 1:shaftNum
        syncSpeedRatio = abs(speedRatio_full(iShaft));
        syncLine = syncSpeedRatio * exciteRad; % Synchronous excitation line for the current shaft
        
        for iRow = 1:size(newEig, 1)
            % Calculate frequency difference
            diffCurve = newEig(iRow, :) - syncLine;
            
            % Find zero-crossing intervals
            crossIdx = find(diffCurve(1:end-1) .* diffCurve(2:end) <= 0);
            
            for k = 1:length(crossIdx)
                idx = crossIdx(k);
                % Extract coordinates before and after the intersection
                x1 = exciteRad(idx);     x2 = exciteRad(idx+1);
                y1 = diffCurve(idx);     y2 = diffCurve(idx+1);
                
                % Avoid division by zero
                if (y2 - y1) == 0
                    xc = x1;
                else
                    % Linear interpolation to find the exact zero-crossing baseline speed
                    xc = x1 - y1 * (x2 - x1) / (y2 - y1);
                end
                
                % Calculate the y-coordinate (natural frequency at the intersection)
                yc = syncSpeedRatio * xc;
                
                % Store the [x, y] coordinates
                criticalPoints = [criticalPoints; xc, yc];
            end
        end
    end
    
    % Remove duplicate coordinate pairs within tolerance and separate the output
    if ~isempty(criticalPoints)
        criticalPoints = uniquetol(criticalPoints, 1e-4, 'ByRows', true);
        criticalSpeeds = criticalPoints(:, 1); % Return only the baseline speeds as requested
    else
        criticalSpeeds = [];
    end
    
    % 5. Plotting
    if options.isPlot
        colors = getDesignPalette();
        figure('Theme', 'Light', 'Name', 'Campbell Diagram', 'Color', colors.Background.Stone.s1);
        hold on; grid on;
        
        % Extract shades from the palette for plotting modes
        accentColors = [colors.Accents.Red.s4; colors.Accents.Blue.s4; ...
                        colors.Accents.Yellow.s4; colors.Extended.Green.s4; ...
                        colors.Extended.Purple.s4; colors.Extended.Teal.s4; ...
                        colors.Extended.Orange.s4; colors.Extended.Olive.s4];
        numColors = size(accentColors, 1);
        
        % Plot natural frequencies
        for iRow = 1:size(newEig, 1)
            plot(exciteRad, newEig(iRow, :), 'Color', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        
        % Plot synchronous excitation lines
        hSyncLines = gobjects(shaftNum, 1);
        for iShaft = 1:shaftNum
            cIdx = mod(iShaft-1, numColors) + 1;
            speedHere = speedRatio_full(iShaft);
            
            hSyncLines(iShaft) = plot(exciteRad, abs(speedHere) * exciteRad, '--', ...
                'Color', accentColors(cIdx,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('%.2fX (\\Omega = %.2f\\omega)', speedHere, speedHere));
        end
        
        % Plot critical speed markers using the 2D coordinate pairs
        if ~isempty(criticalPoints)
            for iCr = 1:size(criticalPoints, 1)
                crX = criticalPoints(iCr, 1);
                crY = criticalPoints(iCr, 2);
                
                % Add the display name only to the first scatter point
                if iCr == 1
                    scatter(crX, crY, 40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors.Accents.Yellow.s3, ...
                        'DisplayName', 'Critical Speeds');
                else
                    scatter(crX, crY, 40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors.Accents.Yellow.s3, ...
                        'HandleVisibility', 'off');
                end
            end
        end
        
        xlabel('Spin Speed \Omega (rad/s)', 'FontSize', 11);
        ylabel('Natural Frequency \omega (rad/s)', 'FontSize', 11);
        title('Campbell Diagram', 'FontSize', 12, 'FontWeight', 'bold');
        
        legend('Location', 'NorthWest');
        
        xlim([min(exciteRad), max(exciteRad)]);
        ylim([0, max(exciteRad)*1.3]);
        box on;
        hold off;
    end
end