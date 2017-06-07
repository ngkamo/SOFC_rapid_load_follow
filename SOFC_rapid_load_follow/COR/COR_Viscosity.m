%% -- SOFC MODEL --
%
function [mu] = COR_Viscosity(T,el)

switch el
    case 'H2'
        mu = 6.162e-6+(1.145e-8)*T;
    case 'H2O'
        mu = 4.567e-6+(2.209e-8)*T; 
    case 'O2'
        mu = 1.668e-5+(3.108e-8)*T;
    case 'CO2'
        mu = 1.273e-5+(2.822e-8)*T;
    case 'CO'
        mu = 1.399e-5+(2.582e-8)*T;
    case 'Air'
        mu = 1.514e-5+(2.793e-8)*T;
    case 'CH4'
        mu = 9.177e-6+(1.755e-8)*T;
    case 'N2'
        mu = 6.162e-6+(1.145e-8)*T;
end

end