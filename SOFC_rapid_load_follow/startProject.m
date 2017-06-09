% ----------------------------------------------
%              Semester Project
% Real-time optimization of a fuell cell system
%         in rapid load following SOFC
%
%     Student: Frederic NGUYEN
% Supervisors: Tafarel DE AVILA FERREIRA
%              Altug BITLISLIOGLU
%   Professor: Dominique BONVIN
% ----------------------------------------------

% Main entry of the program

%%
clear all;
close all;
clc;

addpath('REFORMER','SOFC','COR','HEX','BURN');

%% ---------------------------
%  Initialization
% ----------------------------
global myProblem Ps_el it u_previous

% Number of cells :
N_c = 6;

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);
SOFC_data_plant   = data_SOFC_plant(prjname_SOFC,N_c);
myProblem = data_myProblem(prjname_myProblem);

% -----------
% Parameters
% Stream temperatures:
T_CH4_in  = 200.0; % [C]  methane stream  
T_H2O_in  = 200.0; % [C]  steam stream
T_air_in  = 30;    % [C]  air stream

T_in = [T_CH4_in T_H2O_in T_air_in] + SOFC_data_nominal.cst.K;

% -----------------------------------
% Initial guess for the optimization
% Initial State :
T_r0    = 550; % [C]   fuel reformer

T_el0   = 650; % [C]   electrolyte
T_i0    = 650; % [C]   interconnect 
T_fuel0 = 650; % [C]   fuel channel
T_air0  = 650; % [C]   air channel 

T_b0    = 1000; % [C]  burner

T_h0    = 650; % [C]   heat exchanger - fuel side
T_c0    = 600; % [C]   heat exchanger - air side

T_0 = [T_r0 T_el0 T_i0 T_fuel0 T_air0 T_b0 T_h0 T_c0] + SOFC_data_nominal.cst.K;

% SCTP - CNTP
R = 8.314462175;    % [J/K/mol] 
P = 100e3;          % [Pa]
T_st = 273.15;      % [K]

q_CH4inp      = 0.376129940461854; % [L/min]     methane flow rate
q_AIRcathinp  = 17.2531441172150; % [L/min]     methane flow rate
Iinp          = 19;               % [A]         current

u_0 = [q_CH4inp,q_AIRcathinp,Iinp];

% Initial guess for the optimization problem
u0 = [u_0 T_0];

%% Optimal values
Pel_opt    = [80 90 100];
Ucell_opt  = 0.7;
eff_opt    = [0.4297   0.4260  0.4226];
inputs_opt = [0.3163   0.3589  0.402;
              12.4936 13.8194 15.1185;
              19.0476 21.4286 23.8095];

profile_setpoint = [2 2 2];  % profile setpoint

Ps_el = [];
Uopt_hist = [];
inputs_opt_hist = [];
eff_opt_hist = [];
for i = 1:size(profile_setpoint,2)
    Ps_el(i) = Pel_opt(profile_setpoint(i));
    inputs_opt_hist(:,i) = inputs_opt(:,profile_setpoint(i));
    eff_opt_hist(i) = eff_opt(profile_setpoint(i));
end

%% -------------------------
% Constraints, program loop

% ub = [27.25,272.53/0.21, 30, 600*ones(1,1)+273.15 800*ones(1,4)+273.15 1600+273.15 1474 1474];
% lb = [1.36E-03,0.01/0.21, 0, 450*ones(1,6)+273.15   200+273.15 200+273.15];
ub = [27.25,272.53/0.21, 30, Inf*ones(1,8)];
lb = [1.36E-03,0.01/0.21, 0, -Inf*ones(1,8)];

Aeq = [];
beq = [];

Lair_upper = 10; 
Lair_lower = 3;
FU_upper   = 0.7;
%

kc = (6e+4)*N_c*R*T_st/(8*P*SOFC_data_nominal.cst.F);

A   = [2*Lair_lower,  -0.21,  0,  zeros(1,8);
      -2*Lair_upper,   0.21,  0,  zeros(1,8);
      -FU_upper,          0,  kc, zeros(1,8)];
B   =  [0;0;0];

myProblem.TC.TimeConstantOff = 0;

last_config = [0.3489 16.8194 20.4286 897.5259 1152.4056 1151.7350 1151.9062 1151.4877 1513.9958 978.1512 965.7702]'; % 90 W nominal model
modif = zeros(3,4);

states_plant_hist = [];
power_plant_hist  = [];
volt_plant_hist   = [];
input_plant_hist  = [];
eff_plant_hist    = [];

for i = 1:size(profile_setpoint,2)
    % ----------------------------
    %  RTO layer
    % ----------------------------
    it = i;
    [u_opt] = RTO_layer(u0,modif,T_in,A,B,lb,ub,SOFC_data_nominal,myProblem.OPT.options);

    % ----------------------------
    %  MPC layer
    % ----------------------------
    tic
    [output_plant, input_plant, last_config, P_plant, U_plant, eff_plant, modif] = MPC_layer_grad_2nd(u_opt',modif,Ps_el(it),Ucell_opt,T_in,SOFC_data_nominal,SOFC_data_plant,A,B,lb,ub,last_config);
    toc
    
    states_plant_hist = [ states_plant_hist output_plant];
    input_plant_hist  = [input_plant_hist input_plant];
    power_plant_hist  = [power_plant_hist P_plant];
    volt_plant_hist   = [volt_plant_hist U_plant];
    eff_plant_hist    = [eff_plant_hist eff_plant];
    
    duration_step(i) = size(output_plant,2);
    
    u0 = last_config';
end

%% Optimal values
Popt_hist = repelem(Ps_el,duration_step);
Uopt_hist = Ucell_opt*ones(1,size(volt_plant_hist,2));
eff_opt_hist = repelem(eff_opt_hist,duration_step);
inputs_opt_hist = repelem(inputs_opt_hist,1,duration_step);

%% Post-processing result
set(0,'DefaultFigureWindowStyle','docked')
close all

Ts = 10; % [s] time sample
time_step = [1:1:size(states_plant_hist,2)]*10/60;

figure('Name','temperatures','Color','w')
plot(time_step, states_plant_hist,'LineWidth',1);
legend({'Reformer','electrolyte','interconnect','fuel channel','air channel','burner','HEX fuel','HEX air'},'Interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
ylabel('Temperature [K]','Interpreter','latex')
title('Temperatures in the system','Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','power','Color','w')
plot(time_step, power_plant_hist,time_step,Popt_hist,'--','LineWidth',1);
xlabel('Time [min]','Interpreter','latex')
ylabel('Power $P_{el}$ [W]','Interpreter','latex')
title('Power delivered by the system','Interpreter','latex')
legend({'plant','optimal'},'Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','efficiency','Color','w')
plot(time_step, eff_plant_hist,time_step,eff_opt_hist,'--','LineWidth',1);
xlabel('Time [min]','Interpreter','latex')
ylabel('Efficiency $\eta$ [-]','Interpreter','latex')
title('Efficiency of the the system','Interpreter','latex')
legend({'plant','optimal'},'Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','voltage','Color','w')
plot(time_step, volt_plant_hist,time_step,Uopt_hist,'--','LineWidth',1);
xlabel('Time [min]','Interpreter','latex')
ylabel('Voltage $U_{cell}$ [V]','Interpreter','latex')
title('Voltage of the cell','Interpreter','latex')
legend({'plant','optimal'},'Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','methane','Color','w')
[ts,ys] = stairs(time_step, input_plant_hist(1,:));
plot(ts, ys,time_step,inputs_opt_hist(1,:),'--','Linewidth',1)
xlabel('Time [min]','Interpreter','latex')
ylabel('Methane flow rate $\dot q_{CH4}$ [L/min]','Interpreter','latex')
title('Methane flow rate','Interpreter','latex')
legend({'plant','optimal'},'Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','air','Color','w')
[ts,ys] = stairs(time_step, input_plant_hist(2,:));
plot(ts, ys,time_step,inputs_opt_hist(2,:),'--','Linewidth',1)
xlabel('Time [min]','Interpreter','latex')
ylabel('Air flow rate $\dot q_{air}$ [L/min]','Interpreter','latex')
title('Air flow rate','Interpreter','latex')
legend({'plant','optimal'},'Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on

figure('Name','current','Color','w')
[ts,ys] = stairs(time_step, input_plant_hist(3,:));
plot(ts, ys,time_step,inputs_opt_hist(3,:),'--','Linewidth',1)
xlabel('Time [min]','Interpreter','latex')
ylabel('Current $I$ [A]','Interpreter','latex')
title('Current','Interpreter','latex')
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',12,'TickLabelInterpreter','latex')
grid on