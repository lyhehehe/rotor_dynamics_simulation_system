function generateJacobian(Parameter, calculateOmegaStr, loadMatrix1Str, loadMatrix2Str, processGNStr)

arguments
    Parameter 
    calculateOmegaStr 
    loadMatrix1Str 
    loadMatrix2Str
    processGNStr 
end
%% 

% check the exist of jacobian matrix and create .m
jacobianFileName = 'jacobianJ21J22.m';
if isfile(jacobianFileName)
    delete(jacobianFileName);
end

fJac = fopen(jacobianFileName,'w');


%% 

% write comments line
comments = [];


% write function start
functionStart = [...
"function [J21, J22] = jacobianJ21J22(tn, yn, dyn, Parameter)";...
" "...
];

fprintf(fJac,'%s\n',comments);
fprintf(fJac,'%s\n',functionStart);


%% 

% write status calculation
fprintf(fJac,'%s\n', calculateOmegaStr);


%%

% write matrix loading
fprintf(fJac,'%s\n', loadMatrix1Str);
fprintf(fJac,'%s\n', loadMatrix2Str);


%%

% write G N calculating
fprintf(fJac,'%s\n', processGNStr);

%% 

% write linear part output
linearOutput = [...
" ";...
"% calculate linear part of J21 and J22";...
"J21_lin = -(K - N);";...
"J22_lin = -(C - G);"
];
fprintf(fJac,'%s\n', linearOutput);


%%

% write Herzian nonlinear force part
% Hertzian contact force
if Parameter.ComponentSwitch.hasHertzianForce

    % write codes in dynamicEquation.m
    processHertzianForce = [
" ";...
"% calculate jacobian for Hertzian force";...
"dHerzian_dyn = hertzianForceJacobian(yn, omega, Parameter.Matrix.HerzianParameter, Parameter.Mesh.dofNum);";...
" ";...    
    ]; % write something
    fprintf(fJac,'%s\n', processHertzianForce);
    hertzJacobianJ21 = ' + dHerzian_dyn'; % write plus something
    hertzJacobianJ22 = '';% Herzian force does not contain dyn
else
    hertzJacobianJ21 = '';
    hertzJacobianJ22 = '';
end

%%

% write custom force part
if Parameter.ComponentSwitch.hasCustom
processCustomForce = [...
 " ";...
 "% calculate jacobian for customize force";...
 "[dfCustom_dyn, dfCustom_ddyn] = Parameter.Custom.jacobian(yn, dyn, tn, omega, domega, ddomega, Parameter);";...
 " ";...
];
    fprintf(fJac,'%s\n', processCustomForce);
    customJacobianJ21 = ' + dfCustom_dyn';
    customJacobianJ22 = ' + dfCustom_ddyn';
else
    customJacobianJ21 = '';
    customJacobianJ22 = '';
end

%% 

% write total Jacobian J21 and J22
totalJ21 = {...
 ' ';...
 '% total J21 ';...
 ['J21 = J21_lin', hertzJacobianJ21, customJacobianJ21, ';'];...
 ' ';...
 };
totalJ21 = cell2string(totalJ21);
fprintf(fJac,'%s\n', totalJ21);

totalJ22 = {...
 ' ';...
 '% total J22 ';...
 ['J22 = J22_lin', hertzJacobianJ22, customJacobianJ22, ';'];...
 ' ';...
 };
totalJ22 = cell2string(totalJ22);
fprintf(fJac,'%s\n', totalJ22);


%%

% write function end
functionEnd = [...
"end";...
" "...
];
fprintf(fJac,'%s\n',functionEnd);

%%

% close file

fclose(fJac);
end