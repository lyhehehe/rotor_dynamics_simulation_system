function Response = calculateUnbalanceResponse(Parameter, speedMatrix, options)
% CALCULATEUNBALANCERESPONSE Calculates steady-state unbalance response and optionally plots orbits.
%
% Usage:
%   R = calculateUnbalanceResponse(Parameter)
%       Computes the response using a single default operating condition based on 
%       Parameter.Status.vmax and Parameter.Status.ratio, without plotting.
%
%   R = calculateUnbalanceResponse(Parameter, mySpeedMatrix)
%       Computes responses across multiple custom operating conditions defined in 
%       mySpeedMatrix, without plotting.
%
%   R = calculateUnbalanceResponse(Parameter, [], 'isPlot', true, 'plotNodeID', [3, 5])
%       Computes the default condition and plots the orbit trajectories for Nodes 3 and 5.
%
%   R = calculateUnbalanceResponse(Parameter, mySpeedMatrix, 'isPlot', true, 'plotConditionIdx', 10)
%       Computes custom conditions and plots the orbit for Node 1 (default) at the 10th condition.
%
% Inputs:
%   Parameter   - Struct containing system matrices and mesh info. 
%                 Note: Parameter.Status.vmax should be a scalar (max speed of Shaft 1).
%   speedMatrix - (Optional) A [shaftNum, numSpeeds] double matrix of operating speeds (rad/s).
%                 Each column represents a specific operating condition. 
%                 If empty [], it automatically generates a single-column matrix.
%
% Name-Value Optional Arguments:
%   isPlot           - (Logical) Whether to plot the orbit trajectory. Default is false.
%   plotNodeID       - (Double Array) Array of Node IDs to plot. Default is 1.
%   plotConditionIdx - (Double Array) Array of speed condition indices to plot. Default is 1.
%
% Output:
%   Response - Complex 3D array of size [dofNum, numSpeeds, shaftNum].
%       Dimension 1 (dofNum):    The degree of freedom index of the system.
%       Dimension 2 (numSpeeds): The operating condition index (column of speedMatrix).
%       Dimension 3 (shaftNum):  The EXCITATION SOURCE index. 
%                                Response(:, i, k) contains the vibration of the entire system 
%                                at the i-th condition, excited EXCLUSIVELY by the unbalance forces 
%                                from all shafts running at the EXACT SAME FREQUENCY as the k-th shaft.
%
% Algorithm Details:
%   1. Gyroscopic Matrix Handling:
%      For each operating condition, the total gyroscopic matrix (G_total) is dynamically assembled 
%      by multiplying the base gyroscopic matrix block of each individual shaft by its corresponding 
%      current spin speed. The system dynamic stiffness is assembled using (C - G_total) based on 
%      the right-hand rule convention.
%   2. Excitation Force Handling:
%      - The unbalance input is converted into a base mass-eccentricity phasor (U_base), 
%        where the Y-direction lags the X-direction by 90 degrees (F_y = -j * F_x).
%      - To optimize computation, the algorithm identifies "Unique Speeds" in each condition. 
%        If multiple shafts rotate at the same frequency, their unbalance base vectors are 
%        superimposed (U_combined) before calculating the final force (F = U_combined * omega^2).
%      - This ensures the dynamic stiffness matrix is inverted only once per unique frequency.

    arguments
        Parameter struct
        speedMatrix double = []
        options.isPlot (1,1) logical = false
        options.plotNodeID (1,:) double = 1
        options.plotConditionIdx (1,:) double = 1
    end

    % 1. Speed Matrix Generation & Validation
    shaftNum = Parameter.Shaft.amount;
    
    % Auto-generate speed matrix if not provided (or passed as [])
    if isempty(speedMatrix)
        vmax = Parameter.Status.vmax; 
        speedMatrix = zeros(shaftNum, 1);
        speedMatrix(1, 1) = vmax; % Reference shaft (Shaft 1)
        
        if shaftNum > 1
            speedRatio = Parameter.Status.ratio; % Length should be shaftNum - 1
            for iShaft = 2:shaftNum
                speedMatrix(iShaft, 1) = vmax * speedRatio(iShaft-1);
            end
        end
    end
    
    if size(speedMatrix, 1) ~= shaftNum
        error('The number of rows in speedMatrix must equal the number of shafts.');
    end
    numSpeeds = size(speedMatrix, 2);

    % 2. Extract System Matrices and Mesh Info
    M = Parameter.Matrix.mass;
    C = Parameter.Matrix.damping;
    K = Parameter.Matrix.stiffness;
    G_base = Parameter.Matrix.gyroscopic;
    
    dofNum = Parameter.Mesh.dofNum;
    dofInterval = Parameter.Mesh.dofInterval;
    Node = Parameter.Mesh.Node;
    unbalanceData = Parameter.Matrix.unbalance; % [NodeID, ShaftID, Magnitude, Phase]
    
    % Get shaft DOF range for gyroscopic matrix assembly
    shaftDof = zeros(shaftNum, 2);
    for iShaft = 1:shaftNum
        IShaftNode = Node([Node.onShaftNo] == iShaft & [Node.isBearing] == false);
        startNode = min([IShaftNode.name]);
        endNode = max([IShaftNode.name]);
        shaftDof(iShaft,:) = [dofInterval(startNode,1), dofInterval(endNode,2)];
    end

    % 3. Pre-assemble Unbalance Base Vector (U_base)
    U_base = zeros(dofNum, shaftNum);
    for i = 1:size(unbalanceData, 1)
        nodeID = unbalanceData(i, 1);
        shaftID = unbalanceData(i, 2);
        mag = unbalanceData(i, 3);
        phase = unbalanceData(i, 4); % Phase is strictly in radians
        
        nodeStartDof = dofInterval(nodeID, 1);
        dofX = nodeStartDof;
        dofY = nodeStartDof + 1;
        
        % Right-hand rule convention: Y lags X by 90 degrees
        U_base(dofX, shaftID) = U_base(dofX, shaftID) + mag * exp(1i * phase);
        U_base(dofY, shaftID) = U_base(dofY, shaftID) - 1i * mag * exp(1i * phase);
    end

    % 4. Pre-allocate Response Output Matrix
    Response = zeros(dofNum, numSpeeds, shaftNum);
    speedTol = 1e-5; 

    % 5. Main Calculation Loop
    for iSpeed = 1:numSpeeds
        currentSpeeds = speedMatrix(:, iSpeed);
        
        % Assemble Total Gyroscopic Matrix
        G_total = G_base; 
        for iShaft = 1:shaftNum
            rng = shaftDof(iShaft,1) : shaftDof(iShaft,2);
            G_total(rng, rng) = currentSpeeds(iShaft) * G_total(rng, rng);
        end
        
        activeShafts = find(currentSpeeds > speedTol);
        if isempty(activeShafts)
            continue; 
        end
        
        uniqueSpeeds = unique(currentSpeeds(activeShafts));
        
        for j = 1:length(uniqueSpeeds)
            omega = uniqueSpeeds(j);
            sameSpeedShafts = find(abs(currentSpeeds - omega) < speedTol);
            
            U_combined = sum(U_base(:, sameSpeedShafts), 2);
            F_total = U_combined * (omega^2);
            
            if max(abs(F_total)) > 1e-12
                % Dynamic Stiffness H = -M*w^2 + j*w*(C-G) + K
                H = -(omega^2) * M + 1i * omega * (C - G_total) + K;
                X = H \ F_total;
                
                for k = 1:length(sameSpeedShafts)
                    shaftIdx = sameSpeedShafts(k);
                    Response(:, iSpeed, shaftIdx) = X;
                end
            end
        end
    end
    
    % 6. Post-Processing: Orbit Plotting Integration
    if options.isPlot
        for pNode = options.plotNodeID
            if pNode > size(dofInterval, 1) || pNode < 1
                warning('Node %d does not exist. Skipping plot.', pNode);
                continue;
            end
            for pCond = options.plotConditionIdx
                if pCond > numSpeeds || pCond < 1
                    warning('Condition index %d is out of bounds. Skipping plot.', pCond);
                    continue;
                end
                plotOrbitLocal(Response, Parameter, speedMatrix, pNode, pCond);
            end
        end
    end
end

%% Local Function for Orbit Plotting
function plotOrbitLocal(Response, Parameter, speedMatrix, targetNodeID, speedIndex)
% Extracts physical time-domain responses and plots the total orbit trajectory
    shaftNum = Parameter.Shaft.amount;
    dofInterval = Parameter.Mesh.dofInterval;
    
    dofX = dofInterval(targetNodeID, 1);
    dofY = dofX + 1;
    
    currentSpeeds = speedMatrix(:, speedIndex);
    activeSpeeds = currentSpeeds(currentSpeeds > 1e-5);
    
    if isempty(activeSpeeds)
        warning('Node %d Condition %d: All shafts are at zero speed. No orbit to plot.', targetNodeID, speedIndex);
        return;
    end
    
    % Adaptive time window: enough cycles of the slowest shaft
    minSpeed = min(activeSpeeds);
    maxTime = 10 * (2 * pi / minSpeed); 
    
    % Adaptive sample rate: high enough resolution for the fastest shaft
    maxSpeed = max(activeSpeeds);
    dt = (2 * pi / maxSpeed) / 150; % 150 points per fastest revolution
    t = 0:dt:maxTime;
    
    x_total = zeros(size(t));
    y_total = zeros(size(t));
    
    % Superimpose physical (real) time-domain displacements from all excitation sources
    for k = 1:shaftNum
        omega = currentSpeeds(k);
        if omega > 1e-5
            X_complex = Response(dofX, speedIndex, k);
            Y_complex = Response(dofY, speedIndex, k);
            
            x_total = x_total + real(X_complex * exp(1i * omega * t));
            y_total = y_total + real(Y_complex * exp(1i * omega * t));
        end
    end
    
    % Plot Generation
    figure('Name', sprintf('Orbit of Node %d', targetNodeID));
    plot(x_total, y_total, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2); % Standard MATLAB blue
    hold on;
    
    % Mark the starting coordinate
    plot(x_total(1), y_total(1), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Start (t=0)');
    
    grid on;
    xlabel('X Displacement (m)', 'FontSize', 11);
    ylabel('Y Displacement (m)', 'FontSize', 11);
    
    % Maintained your exact bounding logic
    all_data = [x_total; y_total];
    max_val = max(abs(all_data), [], 'all'); 
    limit_val = max_val * 1.1; 
    
    xlim([-limit_val, limit_val]);
    ylim([-limit_val, limit_val]);
    axis square;
    
    % Convert speeds to RPM for the title string
    rpmStr = sprintf('%.0f ', currentSpeeds * 30/pi);
    title({sprintf('Orbit Trajectory of Node %d', targetNodeID), ...
           sprintf('Condition %d | Speeds: [%s] RPM', speedIndex, strtrim(rpmStr))}, ...
          'FontSize', 12, 'FontWeight', 'bold');
          
    xline(0, 'k--', 'HandleVisibility', 'off');
    yline(0, 'k--', 'HandleVisibility', 'off');
    hold off;
end