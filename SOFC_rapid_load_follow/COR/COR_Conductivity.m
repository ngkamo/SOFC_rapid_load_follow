%% -- SOFC MODEL --
%
function [k] = COR_Conductivity(T,el)


switch el
    case 'H2'
        k = 0.08525+2.964e-4*T;
    case 'H2O'
        k = -0.01450 + 9.782e-5*T;
    case 'O2'
        k = 0.01569 + 5.690e-4*T;
    case 'CO2'
        k = 0.005485 + 6.272e-5*T;
    case 'CO'
        k = 0.01275 + 5.384e-5*T;
    case 'Air'
        k = 0.01329 + 5.539e-5*T;
    case 'CH4'
        k = -0.04446 + 2.093e-4*T;
    case 'N2'
        k = 0.01258 + 5.444e-5*T;
end
