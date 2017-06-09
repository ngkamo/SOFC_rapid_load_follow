% RTO layer with constraint adaptation scheme
function [optimal_var] = RTO_layer(u0,modifiers,T_in,A,B,lb,ub,SOFC_data,opts)
global Ps_el

[optimal_var,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,modifiers,T_in,A,B,lb,ub,SOFC_data,opts);
[~,Ucell_opt,~,~,~,~,~] = fPrimeMyProject(0,optimal_var(4:end),optimal_var(1:3),T_in,SOFC_data);

end