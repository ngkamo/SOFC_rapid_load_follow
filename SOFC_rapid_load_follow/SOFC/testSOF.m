clear all;
close all;
clc;

% addpath('SOFC');

global SOFC_data_nominal myProblem

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,60);
myProblem = data_myProblem(prjname_myProblem);



u = [0.00422798880000000,0.0420000000000000,23];


tspan = [0:3:3000];



inflow  = [0.00754686489821613,0.00636300779089179,0.00171244069544604,0.00251532751999972,2.20584554242723e-07];

Tc_in = 1.0229e+03;
Tref = 856.8358;


myProblem.TC.TimeConstantOff = 0;


x0 = [0.4189    0.4187    0.4188    0.418]*1e+03;

tic
[t_ode, y_ode] = ode15s(@(t, x_SOFC) SOFC_Dynamics(t,x_SOFC,u,inflow,[Tref Tc_in],SOFC_data_nominal), tspan, x0, []);% myProblem.ODE.options
toc


%%
plot(t_ode/60,y_ode-273.15)
xlabel('time [min]','FontSize',14);
ylabel('T [{\circ}C]' ,'FontSize',14);
legend({'Electrolyte','Interconnect','Air','Fuel'},'FontSize',14)

