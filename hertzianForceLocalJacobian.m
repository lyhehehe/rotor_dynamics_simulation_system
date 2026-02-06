function K = hertzianForceLocalJacobian(x, y, thetai, thetao, nb, ri, ro, delta0, kHertz, n)
%HERTZIANFORCELOCALJACOBIAN Calculates the local 2x2 stiffness matrix.

c1 = 2*pi/nb;
c2 = (thetao*ro + thetai*ri) / (ro + ri);

% Initialize stiffness accumulators
k11 = 0; % Sum for xx
k12 = 0; % Sum for xy (symmetric)
k22 = 0; % Sum for yy

for ik = 1:1:nb
    thetak = c1 * (ik-1) + c2; 
    c5 = cos(thetak);
    c6 = sin(thetak);
    deltak = x*c5 + y*c6 - delta0; 
    
    if deltak > 0
        % Derivative of deltak^n is: n * deltak^(n-1)
        stiff_val = n * (deltak^(n-1)); 
        
        % Reuse c5, c6 for matrix elements
        k11 = k11 + stiff_val * c5 * c5;
        k12 = k12 + stiff_val * c5 * c6;
        k22 = k22 + stiff_val * c6 * c6;
    end 
end

% Assemble local matrix
K = kHertz * [k11, k12; k12, k22];

end