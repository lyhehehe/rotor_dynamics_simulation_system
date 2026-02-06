function fHertz = hertzianForce(qn, omega, HerzianParameter, dofNum)

% load parameters
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

omega = [omega, 0]; % add 1 row for adapting the bearing connecting the ground
qn = [qn; 0; 0]; % add 2 rows for adapting the bearing connecting the ground

fHertz2 = zeros(dofNum+2,1); % add 2 rows for adapting the bearing connecting the ground

for iHertz = 1:1:hertzianNum
    x = qn(hertzDof(iHertz,1)) - qn(hertzDof(iHertz,2)); % displacement of the Inner shaft - that of the Outer shaft
    y = qn(hertzDof(iHertz,1)+1) - qn(hertzDof(iHertz,2)+1);
    f = hertzianForceEq(x, y, omega(omegaiNo(iHertz)), omega(omegaoNo(iHertz)), nb(iHertz), ri(iHertz), ro(iHertz), delta0(iHertz), kHertz(iHertz), n(iHertz));
    fHertz2(hertzDof(iHertz,1):hertzDof(iHertz,1)+1) = -f; % the relative displacement is (x_Inner - x_Outer), so the force on the Inner shaft should add a minus "-".
    fHertz2(hertzDof(iHertz,2):hertzDof(iHertz,2)+1) = f; % the force on the outer shaft/ground
end

fHertz = fHertz2(1:end-2);

end