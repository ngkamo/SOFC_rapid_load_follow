function [Entalpy] = energy(T,T_ref,ind)

R = 8.3144621e3;%[J/kmol/K]

% Entalpy = [CH4 H2O H2 CO2 CO]
if ind == 1
Entalpy = (1515+R)*(T-T_ref)+83.47/2*(T^2-T_ref^2)-0.02053/3*(T^3-T_ref^3);
end
if ind == 2
Entalpy = (20750+R)*(T-T_ref)+12.15/2*(T^2-T_ref^2);
end
if ind == 3
Entalpy = (18290+R)*(T-T_ref)+3.719/2*(T^2-T_ref^2);
end
if ind == 4
Entalpy = (23690+R)*(T-T_ref)+31.01/2*(T^2-T_ref^2)-(8.875e-3)/3*(T^3-T_ref^3);
end
if ind == 5
Entalpy = (19660+R)*(T-T_ref)+5.019/2*(T^2-T_ref^2);
end 

Entalpy=Entalpy/1000;%[J/mol]

end