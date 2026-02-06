
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
G(1:84, 1:84) = domega(1)*G(1:84, 1:84);
N(1:84, 1:84) = ddomega(1)*N(1:84, 1:84);
G(85:152, 85:152) = domega(2)*G(85:152, 85:152);
N(85:152, 85:152) = ddomega(2)*N(85:152, 85:152);
 
% calculate linear part of J21 and J22
J21_lin = -(K - N);
J22_lin = -(C - G);
 
% calculate jacobian for Hertzian force
dHerzian_dyn = hertzianForceJacobian(yn, omega, Parameter.Matrix.HerzianParameter, Parameter.Mesh.dofNum);
 
 
% calculate jacobian for customize force
[dfCustom_dyn, dfCustom_ddyn] = Parameter.Custom.jacobian(yn, dyn, tn, omega, domega, ddomega, Parameter);
 
 
% total J21 
J21 = J21_lin + dHerzian_dyn + dfCustom_dyn;
 
 
% total J22 
J22 = J22_lin + dfCustom_ddyn;
 
end
 
