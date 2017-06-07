%% -- SOFC MODEL --
%
%  [s] = SOFC_Entropy(T_fuel,T_air)
%  [s] = [sH2an, sH2Oan, sO2cath, sN2cath];
%
function [s] = COR_Entropy(T_fuel,T_air,x_out)

R = 8.314462175;

if x_out(1) <=0
   x_out(1) = 1e-200; 
end

sH2Oan = 188.83 + 143.05*log(T_fuel/298.15) - 232.16*(T_fuel^.25 - 298.15^.25) + 16.5502*(T_fuel^.5 - 298.15^.5) - .036989*(T_fuel - 298.15) - R*log(x_out(2));
sH2an = 130.59 + 56.505*log(T_fuel/298.15) + (22222.6/.75)*(T_fuel^(-.75) - 298.15^(-.75)) - 116500*(T_fuel^(-1) - 298.15^(-1)) + 373800*(T_fuel^(-1.5) - 298.15^(-1.5)) - R*log(x_out(1));
sO2cath = 205.14 + 37.432*log(T_air/298.15) + 2.0102e-5/1.5*(T_air^1.5 - 298.15^1.5) + (178570/1.5)*(T_air^(-1.5)-298.15^(-1.5)) - 1184400*(T_air^(-2) - 298.15^(-2)) - R*log(x_out(3));
sN2cath = 191.61 + 39.06*log(T_air/298.15) + 3.4186e5*(T_air^(-1.5) - 298.15^(-1.5)) - 5363.5e3*(T_air^(-2) - 298.15^(-2)) + 2.73467e8*(T_air^(-3) - 298.15^(-3)) - R*log(x_out(4));

s = [sH2an, sH2Oan, sO2cath, sN2cath];

end

