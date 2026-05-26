
function [J21, J22] = jacobianJ21J22(tn, yn, dyn, Parameter)
 
 
% calculate phase, speed and acceleration
if tn <= 10
    ddomega = [20  26];
    domega  = [20  26] * tn;
    omega   = 0.5 * [20  26] * tn^2;
elseif tn <= 10
    ddomega = [0  0];
    domega  = [200  260];
    omega   = [1000  1300] + [200  260] * (tn - 10 );
elseif tn <= 20
    ddomega = -[20  26];
    domega  = [200  260] - [20  26] * (tn - 10 );
    omega   = [1000  1300] + [200  260]*(tn - 10 ) - 0.5*[20  26]*(tn - 10 )^2;
else
    ddomega = [0  0];
    domega  = [0  0];
    omega   = [2000  2000] + [0  0]*(tn - 20 );
end
 
 
% load matrix
G = Parameter.Matrix.gyroscopic;
N = Parameter.Matrix.matrixN;
K = Parameter.Matrix.stiffness;
C = Parameter.Matrix.damping;
G(1:68, 1:68) = domega(1)*G(1:68, 1:68);
N(1:68, 1:68) = ddomega(1)*N(1:68, 1:68);
G(69:116, 69:116) = domega(2)*G(69:116, 69:116);
N(69:116, 69:116) = ddomega(2)*N(69:116, 69:116);
 
% calculate linear part of J21 and J22
J21_lin = -(K - N);
J22_lin = -(C - G);
 
% total J21 
J21 = J21_lin;
 
 
% total J22 
J22 = J22_lin;
 
end
 
