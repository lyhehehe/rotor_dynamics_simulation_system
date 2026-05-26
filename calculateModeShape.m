%% calculateModeShape - Calculate and visualize geometric mode shapes for a multi-shaft rotor system
%
% This function extracts synchronous mode shape vectors at specified shaft speeds
% and optionally renders the deformed rotor geometry with shaft cross-sections.
%
%% Syntax
%  [ModeShapes, ZCoords] = calculateModeShape(Parameter, exciteRad)
%  [ModeShapes, ZCoords] = calculateModeShape(Parameter, exciteRad, Name=Value)
%
%% Description
% |calculateModeShape| assembles the state-space form, computes eigenvalues and
% eigenvectors, and extracts the translational mode shape at each requested speed.
% The function:
% * Builds and solves the first-order eigenvalue problem at each speed in |exciteRad|
% * Selects the mode whose natural frequency is closest to the current spin speed
% * Accounts for shaft offsets introduced by intermediate bearings when computing global Z-coordinates
% * Renders deformed shaft geometry with spline-smoothed cross-sections if |isPlot| is true
%
%% Input Arguments
% * |Parameter| - System configuration structure containing:
%   * |Matrix|: Global matrices (mass, stiffness, damping, gyroscopic)
%   * |Mesh|: Mesh data including |dofNum|, |dofInterval|, |nodeDistance|, |ShaftElement|
%   * |Status|: Speed ratio |ratio| for multi-shaft gyroscopic scaling
%   * |Shaft|: Shaft properties including |amount|
% * |exciteRad| - Target spin speeds [rad/s] at which to extract mode shapes [1×N double]
%
%% Name-Value Arguments
% * |isPlot|          - Render deformed geometry plots (default: false) [logical]
% * |direction|       - Translational direction to extract: 'X' or 'Y' (default: 'X') [char]
% * |modeSelect|      - Mode selection strategy: 'sync' — closest to excitation frequency (default: 'sync') [char]
% * |scaleFactor|     - Deformation scale: 'auto' or a numeric multiplier (default: 'auto') [char|double]
% * |isUseGyroMatrix| - Include gyroscopic matrix in eigenvalue problem (default: true) [logical]
%
%% Output Arguments
% * |ModeShapes| - Mode shape displacement vectors, one cell per shaft per speed [shaftNum×N cell]
% * |ZCoords|    - Global Z-coordinates of shaft nodes accounting for intermediate bearing offsets [shaftNum×1 cell]
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function [ModeShapes, ZCoords] = calculateModeShape(Parameter, exciteRad, options)
    arguments
        Parameter struct
        exciteRad (1,:) double
        options.isPlot (1,1) logical = false
        options.direction (1,:) char {mustBeMember(options.direction, {'X', 'Y'})} = 'X'
        options.modeSelect (1,:) char {mustBeMember(options.modeSelect, {'sync'})} = 'sync'
        options.scaleFactor = 'auto'
        options.isUseGyroMatrix (1,1) logical = true % added option
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
    
    % Map direction to local DOF offset (1 for X, 2 for Y)
    if strcmp(options.direction, 'X')
        dirOffset = 0;
    else
        dirOffset = 1;
    end
    
    Node = Parameter.Mesh.Node;
    dofInterval = Parameter.Mesh.dofInterval;
    shaftDof = zeros(shaftNum, 2);
    for iShaft = 1:shaftNum
        IShaftNode = Node([Node.onShaftNo] == iShaft & [Node.isBearing] == false);
        startNode = min([IShaftNode.name]);
        endNode = max([IShaftNode.name]);
        shaftDof(iShaft,:) = [dofInterval(startNode,1), dofInterval(endNode,2)];
    end
    
    % Precompute matrix divisions outside the loop
    invM_K = M \ K;
    invM_C = M \ C;
    
    % Arrays to store extracted vectors
    rawEigVectors = zeros(eigNum, exciteRadNum);
    
    % 2. Main Calculation Loop
    if options.isUseGyroMatrix
        % --- With gyroscopic effect: reassemble matrix and solve eigenvalue at each iteration ---
        for iRad = 1:exciteRadNum
            iBasicSpeed = exciteRad(iRad);
            G = G_base; 
            
            for iShaft = 1:shaftNum
                if iShaft == 1
                    shaftSpeed = iBasicSpeed;
                else
                    shaftSpeed = iBasicSpeed * speedRatio(iShaft-1);
                end
                rng = shaftDof(iShaft,1) : shaftDof(iShaft,2);
                G(rng, rng) = shaftSpeed * G(rng, rng);
            end
            
            invM_G = M \ G;
            A = [-(invM_C + invM_G), -invM_K; ...
                 eye(dofNum),        zeros(dofNum, dofNum)];
            
            [V, D] = eig(full(A));
            eigenvalues = diag(D);
            
            % Mode Selection Strategy ('sync': closest to excitation frequency)
            targetFreq = iBasicSpeed; 
            [~, targetIdx] = min(abs(abs(imag(eigenvalues)) - targetFreq));
            
            % Extract complex mode shape vector
            trans = V(:, targetIdx);
            rawEigVectors(:, iRad) = real(trans) + imag(trans); 
        end
        
    else
        % --- Without gyroscopic effect: state matrix is constant, solve eigenvalue only once ---
        A = [-invM_C,     -invM_K; ...
             eye(dofNum), zeros(dofNum, dofNum)];
             
        [V, D] = eig(full(A));
        eigenvalues = diag(D);
        
        % Still loop over each speed to match the closest frequency mode shape from the constant modal pool
        for iRad = 1:exciteRadNum
            targetFreq = exciteRad(iRad); 
            [~, targetIdx] = min(abs(abs(imag(eigenvalues)) - targetFreq));
            
            trans = V(:, targetIdx);
            rawEigVectors(:, iRad) = real(trans) + imag(trans); 
        end
    end
    
    % 3. Calculate Global Node Positions (Handling Intermediate Bearings)
    nodeDistance = Parameter.Mesh.nodeDistance;
    offsetPosition = zeros(shaftNum, 1);
    
    if isfield(Parameter, 'IntermediateBearing')
        InterBearing = Parameter.IntermediateBearing;
        for iShaft = 2:shaftNum
            for iInterBearing = 1:InterBearing.amount
                if InterBearing.betweenShaftNo(iInterBearing, 2) == iShaft
                    basicShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 1);
                    laterShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 2);
                    basicShaftNo = InterBearing.betweenShaftNo(iInterBearing, 1);
                    offsetPosition(iShaft) = basicShaftD - laterShaftD + offsetPosition(basicShaftNo);
                end 
            end 
        end 
    end 
    
    ZCoords = cell(shaftNum, 1);
    for iShaft = 1:shaftNum
        ZCoords{iShaft} = nodeDistance{iShaft} + offsetPosition(iShaft);
    end
    
    % 4. Extract Translation Displacement for Specific Direction
    targetDof = zeros(length(Node), 1);
    counter = 1;
    for iNode = 1:length(Node)
        if ~Node(iNode).isBearing
            targetDof(counter) = dofInterval(iNode, 1) + dirOffset;
            counter = counter + 1;
        end
    end
    targetDof(targetDof == 0) = [];
    eigVectorPart = rawEigVectors(targetDof, :);
    
    % Divide eigen vector according to shaft No.
    ModeShapes = cell(shaftNum, exciteRadNum);
    for iRad = 1:exciteRadNum
        counter = 1;
        for iShaft = 1:shaftNum
            nodeNumThisShaft = length(nodeDistance{iShaft});
            ModeShapes{iShaft, iRad} = eigVectorPart(counter:counter+nodeNumThisShaft-1, iRad);
            counter = counter + nodeNumThisShaft;
        end
    end
    
    % 5. Deformed Geometry Render Engine
    if options.isPlot
        colors = getDesignPalette();
        accentColors = [colors.Accents.Blue.s4; colors.Accents.Red.s4; ...
                        colors.Accents.Yellow.s4; colors.Extended.Green.s4; ...
                        colors.Extended.Purple.s4; colors.Extended.Teal.s4];
        numColors = size(accentColors, 1);
        
        ShaftElement = Parameter.Mesh.ShaftElement;
        
        % Find maximum outer radius across all shafts for auto-scaling
        maxRo = 0;
        for k = 1:length(ShaftElement)
            if ShaftElement(k).outerRadius > maxRo
                maxRo = ShaftElement(k).outerRadius;
            end
        end
        if maxRo == 0, maxRo = 0.05; end % Fallback
        
        for iRad = 1:exciteRadNum
            figName = sprintf('Mode Shape at \\Omega = %.1f rad/s', exciteRad(iRad));
            hFig = figure('Theme', 'Light', 'Name', figName, 'Color', colors.Background.Stone.s1, ...
                          'Units', 'normalized', 'Position', [0.15, 0.2, 0.7, 0.5]);
            ax = axes('Parent', hFig);
            hold(ax, 'on'); grid(ax, 'on');
            
            % Determine Scale Factor
            maxDisp = 0;
            for iShaft = 1:shaftNum
                maxDisp = max(maxDisp, max(abs(ModeShapes{iShaft, iRad})));
            end
            if ischar(options.scaleFactor) && strcmp(options.scaleFactor, 'auto')
                if maxDisp > 0
                    % Scale the maximum deflection to be exactly 1.5 times the maximum radius
                    scaleFactor = (1.5 * maxRo) / maxDisp;
                else
                    scaleFactor = 1;
                end
            else
                scaleFactor = options.scaleFactor;
            end
            
            hLegendLines = gobjects(shaftNum, 1);
            
            nodeNumBefore = 0;
            for iShaft = 1:shaftNum
                cIdx = mod(iShaft-1, numColors) + 1;
                shaftColor = accentColors(cIdx, :);
                
                zNodes = ZCoords{iShaft};
                uNodes = ModeShapes{iShaft, iRad} * scaleFactor;
                
                % Spline interpolation for smooth center line
                zSmooth = linspace(min(zNodes), max(zNodes), 300);
                uSmooth = spline(zNodes, uNodes, zSmooth);
                
                currentShaftElems = ShaftElement([ShaftElement.ShaftNo] == iShaft);
                
                % Plot each element as a deformed patch
                if iShaft > 1
                    nodeNumBefore = nodeNumBefore + length(ZCoords{iShaft-1});
                end
                
                for k = 1:length(currentShaftElems)
                    elem = currentShaftElems(k);
                    Ro = elem.outerRadius;
                    Ri = elem.innerRadius;
                    
                    % Get element start and end global Z coordinates
                    zStart = zNodes(elem.NodeNo(1)-nodeNumBefore);
                    zEnd = zNodes(elem.NodeNo(2)-nodeNumBefore);
                    
                    % Filter smooth points belonging to this element
                    idx = (zSmooth >= zStart) & (zSmooth <= zEnd);
                    zElem = zSmooth(idx);
                    uElem = uSmooth(idx);
                    
                    % Ensure exact endpoints are included to prevent gaps
                    if zElem(1) > zStart
                        zElem = [zStart, zElem]; uElem = [spline(zNodes, uNodes, zStart), uElem];
                    end
                    if zElem(end) < zEnd
                        zElem = [zElem, zEnd]; uElem = [uElem, spline(zNodes, uNodes, zEnd)];
                    end
                    
                    % --- Draw Deformed Shaft Walls ---
                    if Ri > 0
                        % Hollow Shaft: Draw Top Wall and Bottom Wall separately
                        topWall_Y = uElem + Ro;
                        topWall_inner_Y = uElem + Ri;
                        patch(ax, [zElem, fliplr(zElem)], [topWall_Y, fliplr(topWall_inner_Y)], ...
                              shaftColor, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
                          
                        botWall_Y = uElem - Ro;
                        botWall_inner_Y = uElem - Ri;
                        patch(ax, [zElem, fliplr(zElem)], [botWall_inner_Y, fliplr(botWall_Y)], ...
                              shaftColor, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
                    else
                        % Solid Shaft: Draw entire section
                        topWall_Y = uElem + Ro;
                        botWall_Y = uElem - Ro;
                        patch(ax, [zElem, fliplr(zElem)], [topWall_Y, fliplr(botWall_Y)], ...
                              shaftColor, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
                    end
                end
                
                % Plot Center Line
                plot(ax, zSmooth, uSmooth, '-.', 'Color', shaftColor * 0.5, 'LineWidth', 1.2, 'HandleVisibility', 'off');
                
                % Plot Nodes
                scatter(ax, zNodes, uNodes, 25, 'filled', 'MarkerFaceColor', colors.Background.Stone.s1, ...
                        'MarkerEdgeColor', shaftColor, 'LineWidth', 1, 'HandleVisibility', 'off');
                
                % Create dummy line for Legend
                hLegendLines(iShaft) = plot(ax, NaN, NaN, '-', 'Color', shaftColor, 'LineWidth', 4, ...
                                            'DisplayName', sprintf('Shaft %d', iShaft));
            end
            
            % Beautify plot with black text colors
            xlabel(ax, 'Global Axis Position Z (m)', 'FontName', 'Helvetica', 'FontSize', 11, 'Color', 'k');
            ylabel(ax, 'Deformation Amplitude (Scaled)', 'FontName', 'Helvetica', 'FontSize', 11, 'Color', 'k');
            title(ax, figName, 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
            
            % Layout adjustments
            globalMinZ = min(cellfun(@min, ZCoords));
            globalMaxZ = max(cellfun(@max, ZCoords));
            xlim(ax, [globalMinZ - 0.1, globalMaxZ + 0.1]);
            
            yLimMax = maxRo * 3; % Provide enough breathing room around the geometries
            ylim(ax, [-yLimMax, yLimMax]);
            
            % Equal data aspect ratio to maintain true physical geometry proportions
            axis(ax, 'equal');
            
            % Default MATLAB legend formatting
            legend(ax, hLegendLines, 'Location', 'best');
            
            set(ax, 'Box', 'on', 'TickDir', 'in', 'LineWidth', 0.8, 'Color', 'w');
        end
    end
end