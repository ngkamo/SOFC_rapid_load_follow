%% -- SOFC MODEL --
%
function [rho] = COR_Density(T,p,el)

switch el
    case 'H2'
        Mm = 2.01588;
    case 'H2O'
        Mm = 18.01528;
    case 'O2'
        Mm = 31.99880;
    case 'CO2'
        Mm = 44.0095;
    case 'CO'
        Mm = 28.0101;
    case 'Air'
        Mm = 28.97;
    case 'CH4'
        Mm = 16.0425;
    case 'N2'
        Mm = 28.01340;
end

R = 8.314462;
Rm = R/Mm*1000;

rho = p/(Rm*T);

end