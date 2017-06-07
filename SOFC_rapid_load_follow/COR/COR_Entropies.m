%% -- SOFC MODEL --
%
function [Ds] = COR_Entropies(T_ref,T,el)

switch el
    case 'H2'
        Ds = 18290*log(T/T_ref)+3.719*(T-T_ref);
    case 'H2O'
        Ds = 20750*log(T/T_ref)+12.15*(T-T_ref); 
    case 'O2'
        Ds = 22010*log(T/T_ref)+4.936*(T-T_ref);
    case 'CO2'
        Ds = 23690*log(T/T_ref)+31.01*(T-T_ref)-(8.875e-3)/2*(T^2-T_ref^2);
    case 'CO'
        Ds = 19660*log(T/T_ref)+5.019*(T-T_ref);
    case 'Air'
        Ds = 19480*log(T/T_ref)+4.936*(T-T_ref);
    case 'CH4'
        Ds = 1515*log(T/T_ref)+83.47*(T-T_ref)-0.02053/2*(T^2-T_ref^2);
    case 'N2'
        Ds = 19100*log(T/T_ref)+5.126*(T-T_ref);
end

Ds=Ds/1000;

end