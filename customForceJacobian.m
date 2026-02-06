function [dfCustomd_dyn, dfCustomd_ddyn] = customForceJacobian(yn, dyn, tn, omega, domega, ddomega, Parameter)

dof_num = Parameter.Mesh.dofNum;

dfCustomd_dyn = sparse(dof_num, dof_num);
dfCustomd_ddyn = dfCustomd_dyn;

end