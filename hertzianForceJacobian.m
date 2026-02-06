function JHertz = hertzianForceJacobian(qn, omega, HerzianParameter, dofNum)
%HERTZIANFORCEJACOBIAN Calculates the Jacobian matrix (tangent stiffness) of Hertzian contact forces.

% Load parameters
nb = HerzianParameter.nb;
ri = HerzianParameter.ri;
ro = HerzianParameter.ro;
delta0 = HerzianParameter.delta0;
kHertz = HerzianParameter.kHertz;
n = HerzianParameter.n;
omegaiNo = HerzianParameter.omegaiNo;
omegaoNo = HerzianParameter.omegaoNo;
hertzDof = HerzianParameter.hertzDof;
hertzianNum = HerzianParameter.hertzianNum;

% Pre-processing input vectors for ground adaptation
omega = [omega, 0]; 
qn = [qn; 0; 0]; 

% Pre-allocate arrays for sparse matrix triplets
% Each bearing contributes to a 4x4 block (16 elements max)
maxNz = hertzianNum * 16;
rows = zeros(maxNz, 1);
cols = zeros(maxNz, 1);
vals = zeros(maxNz, 1);
count = 1;

for iHertz = 1:1:hertzianNum
    % Extract indices
    idxInX = hertzDof(iHertz,1);
    idxInY = hertzDof(iHertz,1)+1;
    idxOutX = hertzDof(iHertz,2);
    idxOutY = hertzDof(iHertz,2)+1;
    
    % Calculate relative displacement (Inner - Outer)
    x = qn(idxInX) - qn(idxOutX);
    y = qn(idxInY) - qn(idxOutY);
    
    % Call local Jacobian kernel
    K_loc = hertzianForceLocalJacobian(x, y, omega(omegaiNo(iHertz)), omega(omegaoNo(iHertz)), ...
        nb(iHertz), ri(iHertz), ro(iHertz), delta0(iHertz), kHertz(iHertz), n(iHertz));
    
    % Assemble 4x4 block into triplets
    % The structure is: [ -K_loc,  K_loc; 
    %                      K_loc, -K_loc ]
    
    % Define Global Indices for this block
    gIdx = [idxInX, idxInY, idxOutX, idxOutY];
    
    % Iterate through the 2x2 local matrix to fill the 4x4 global block
    for r = 1:2
        for c = 1:2
            val = K_loc(r,c);   
            % Block (1,1): Inner-Inner (Partial Force_In / Partial Disp_In) -> -K_loc
            rows(count) = gIdx(r); cols(count) = gIdx(c); vals(count) = -val; count = count + 1;
            
            % Block (1,2): Inner-Outer (Partial Force_In / Partial Disp_Out) -> +K_loc
            rows(count) = gIdx(r); cols(count) = gIdx(c+2); vals(count) = val; count = count + 1;
            
            % Block (2,1): Outer-Inner (Partial Force_Out / Partial Disp_In) -> +K_loc
            rows(count) = gIdx(r+2); cols(count) = gIdx(c); vals(count) = val; count = count + 1;
            
            % Block (2,2): Outer-Outer (Partial Force_Out / Partial Disp_Out) -> -K_loc
            rows(count) = gIdx(r+2); cols(count) = gIdx(c+2); vals(count) = -val; count = count + 1;
        end
    end
end

% Filter out valid entries (remove unused pre-allocated slots)
validIdx = 1:(count-1);
rows = rows(validIdx);
cols = cols(validIdx);
vals = vals(validIdx);

% Filter out Ground DOFs (indices > dofNum)
mask = (rows <= dofNum) & (cols <= dofNum);
JHertz = sparse(rows(mask), cols(mask), vals(mask), dofNum, dofNum);

end