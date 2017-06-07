function [myf,myc,myceq] = computeall(u0,modif,T_in,SOFC_data)
global  Ps_el it u_previous

if isempty(modif) == 1
   modif = zeros(3,4);
end

if isempty(it) == 1
   it = 1;
end

if isempty(u_previous) == 1
    u_previous = [0 0 0];
end

u = u0(1:3);
x_steady = u0(4:end);

[dx_global,U_c,P_el,~,~,~,eta_sys] = fPrimeMyProject(0,x_steady,u,T_in,SOFC_data);


myc = [];
myc(1) = 0.7 - U_c - modif(2,4) - modif(2,1:3)*(u-u_previous)';              % Cell voltage

myceq = [P_el-Ps_el(it) + modif(3,4) + modif(3,1:3)*(u-u_previous)';       % Power
         dx_global];

myf = - eta_sys - modif(1,1:3)*u';

end