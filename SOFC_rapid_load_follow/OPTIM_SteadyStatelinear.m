function [tempsteady_SOFC] = OPTIM_SteadyStatelinear(u,x0,T_in,SOFC_data)
opts = optimset('Diagnostics','off', 'Display','off');
[tempsteady_SOFC,FVAL,EXIT] = fsolve(@(x) fPrimeMyProject(0,x,u,x_imposed,T_in,SOFC_data),x0,opts);

end