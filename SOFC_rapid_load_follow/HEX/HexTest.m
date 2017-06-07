%% Hex test
clear all
clc

close all
% global  Th_in Tc_in inflowh inflowc
global SOFC_data_nominal myProblem
addpath('COR')
prjname_SOFC       = 'Solid Oxide Fuel Cell';
N_c = 6;
SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);

prjname_myProblem  = 'Real Time Optimization';
myProblem = data_myProblem(prjname_myProblem);

% Initial conditions

% x0 = [700;600;500];
% x0 = [1021.28183017018,903.368241151310];
t_start = 0;
t_final = 5000;

% options = optimset('Algorithm','trust-region-reflective',...
%                    'LargeScale','on',...
%                    'Diagnostics','on',...
%                    'Display','off' ,...
%                    'FunValCheck','on',...
%                    'MaxFunEvals',10000,...
%                    'MaxIter',10000,...
%                    'TolFun',1e-12,...
%                    'TolX',1e-12);

Th_in = 1.503511797052666e+03; % temperature of the burnt fuel from the burner
Tc_in = 30+273.15; % temperature of the air
% inflowh =  [0 0 0.0021 0.0003 0 0.0018 0.0113]; % molar flow of the stream coming from the burner
% inflowc = [0.0030 0.0113];                      % molar flow of the air
% Equilibrium

u = [0.000157727181360379/10,0.00100407731574690/10,13.1064285367178];
inflowc = [u(2) u(2)*SOFC_data_nominal.CH.r_N2];

inflow = [0,0.0173347540800000,0,0,0.00422798880000000,0.0335440224000000,0.157999986886193,0.0420000000000000];
% inflow = [inflowh inflowc];
% T_HEX      = [x_BUR TO2in];
T_HEX      = [Th_in Tc_in];

% [dxH, outflow] = HEX_Dynamics(xH,inflow,T_HEX,SOFC_data_nominal)
% tic
% xSS = fsolve(@(x) HEX_Dynamics(0,x,SOFC_data_nominal),x0, options)
% [xSS,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) HEX_Dynamics(0,x,inflow,T_HEX,SOFC_data_nominal),x0, options)
% toc
% xSS-273.15




% Dynamic simulation
% x0 = [1200;300;600];
x0 = [700 30]+273.15;
tspan = t_start:1:t_final;

disp(' ');
disp(' Solving ODE 45 ');
disp(' ');

tic
[t, y] = ode15s(@(t, x) HEX_Dynamics(t, x,inflow,T_HEX,SOFC_data_nominal), tspan, x0);
toc                                  
y(end,:)
% % Plotting

% plot(t/60,y(:,1),t/60,y(:,2),t/60,y(:,3));
plot(t/60,y(:,1),t/60,y(:,2));
% legend('T_hout','T_cout','T_{met}');
legend('T_hout','T_cout');
grid on;
title('Evolution of temperature through time')
xlabel('Time [min]')
ylabel('Temperature [K]')