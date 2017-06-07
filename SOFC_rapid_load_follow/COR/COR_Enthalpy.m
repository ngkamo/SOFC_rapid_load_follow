%% -- SOFC MODEL --
%
%  [h] = SOFC_Enthalpy(T_fuel,T_air)
%  [h] = [hH2an hH2Oan hO2cath hN2cath];
%
function [h] = COR_Enthalpy(T_fuel,T_air)

hH2Oan  = -241826 + 143.05*(T_fuel - 298.15) - 46.432*(T_fuel^1.25 - 298.15^1.25) + (8.2751/1.5)*(T_fuel^1.5 - 298.15^1.5) - .0184945*(T_fuel^2 - 298.15^2);
hH2an   = 56.505*(T_fuel - 298.15) - 88890.4*(T_fuel^.25 - 298.15^.25) + 116500*log(T_fuel/298.15) + 1121400*(T_fuel^(-.5) - 298.15^(-.5));
hO2cath = 37.432*(T_air - 298.15) + 8.0408e-6*(T_air^2.5 - 298.15^2.5) + 357140*(T_air^(-.5) - 298.15^(-.5)) - 2368800*(T_air^(-1) - 298.15^(-1));
hN2cath = 39.06*(T_air - 298.15) + 1025580*(T_air^(-.5)-298.15^(-.5)) - 10727e3*(T_air^(-1)-298.15^(-1)) + 4.102e8*(T_air^(-2) - 298.15^(-2));

h = [hH2an hH2Oan hO2cath hN2cath];

end

