%% -- SOFC MODEL --
%
function [Cp] = COR_CP(T,el)

R = 8.3144621e3;%[J/kmol/K]

switch el
    case 'H2'
        Cp = (18290+R)+3.719*T;
    case 'H2O'
        Cp = (20750+R)+12.15*T; 
    case 'O2'
        Cp = (22010+R)+4.936*T;
    case 'CO2'
        Cp = (23690+R)+31.01*T-8.875e-3*T^2;
    case 'CO'
        Cp = (19660+R)+5.019*T;
    case 'Air'
        Cp = (19480+R)+4.936*T;
    case 'CH4'
        Cp = (1515+R)+83.47*T-0.02053*T^2;
    case 'N2'
        Cp = (19100+R)+5.126*T;
end

Cp=Cp/1000;

% t=T/1000;
% 
% %CH4 298->1300K
% 
% A = -0.703029;
% B = 108.4773;
% C = -42.52157;
% D = 5.862788;
% E = 0.678565;
% F = -76.84376;
% G = 158.7163;
% H = -74.87310;
% 
% Cp = A + B*t + C*t.^2 + D*t.^3 + E./t.^2;
% Dh = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
% S  = A*ln(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G;

end