%% -- SOFC MODEL --
%
function [Dh] = COR_Energies(T_ref,T,el)

R = 8.3144621e3;%[J/kmol/K]

switch el
    case 'H2'
        Dh = (18290+R)*(T-T_ref)+3.719/2*(T^2-T_ref^2);
    case 'H2O'
        Dh = (20750+R)*(T-T_ref)+12.15/2*(T^2-T_ref^2); 
    case 'O2'
        Dh = (22010+R)*(T-T_ref)+4.936/2*(T^2-T_ref^2);
    case 'CO2'
        Dh = (23690+R)*(T-T_ref)+31.01/2*(T^2-T_ref^2)-(8.875e-3)/3*(T^3-T_ref^3);
    case 'CO'
        Dh = (19660+R)*(T-T_ref)+5.019/2*(T^2-T_ref^2);
    case 'Air'
        Dh = (19480+R)*(T-T_ref)+4.936/2*(T^2-T_ref^2);
    case 'CH4'
        Dh = (1515+R)*(T-T_ref)+83.47/2*(T^2-T_ref^2)-0.02053/3*(T^3-T_ref^3);
    case 'N2'
        Dh = (19100+R)*(T-T_ref)+5.126/2*(T^2-T_ref^2);
end

Dh=Dh/1000;

end