%% Burn test
clear all
close all
clc
warning on
global  Tf_in Ta_in inflowf inflowa SOFC_data_nominal myProblem
addpath('COR')
prjname_SOFC       = 'Solid Oxide Fuel Cell';
N_c = 6;
SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);
prjname_myProblem  = 'Real Time Optimization';
myProblem = data_myProblem(prjname_myProblem);

% Initial conditions

x0 = 1000;
t_start = 0;
t_final = 1800;

options = optimset('Algorithm','trust-region-reflective',...
                   'LargeScale','on',...
                   'Diagnostics','on',...
                   'Display','off' ,...
                   'FunValCheck','on',...
                   'MaxFunEvals',10000,...
                   'MaxIter',10000,...
                   'TolFun',1e-12,...
                   'TolX',1e-12);

Tf_in = 1.0787e+03;
Ta_in = 1.0790e+03;
inflowf = [0.00724533176340312,0.0100894223165969,0,0.00251532687543774,0.00171266192456226,0.0384243517194204,0.157999986886193]/5;
Tref  = 560+273.15;
% Equilibrium

%x0 = fsolve(@(x) BURN_Dynamics(0,x),x0, options);

% Dynamic simulation

tspan = t_start:0.01:t_final;

disp(' ');
disp(' Solving ODE ');
disp(' ');

tic
[t, y] = ode15s(@(t, x) BURN_Dynamics(t, x, Tf_in, Ta_in, inflowf,Tref,SOFC_data_nominal), tspan, x0, options);
toc                                          
    
% Plotting

plot(t/60,y(:,1)-273.15);
legend('T_b');
grid on;
title('Evolution of temperature through time')
xlabel('Time [min]')
ylabel('Temperature [K]')