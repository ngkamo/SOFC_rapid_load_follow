clear all;
close all;
clc;
warning on
addpath('REFORMER','SOFC','COR','HEX','BURN');

%%
global myProblem Ps_el

% Number of cells :
N_c = 6;

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_plant(prjname_SOFC,N_c);
myProblem = data_myProblem(prjname_myProblem);

%%   Parameters
% Stream temperatures:
T_CH4_in  = 200.0; % [C]  methane stream  
T_H2O_in  = 200.0; % [C]  steam stream
T_air_in  = 30;    % [C]  air stream

T_in = [T_CH4_in T_H2O_in T_air_in] + SOFC_data_nominal.cst.K;

%% Initial guess for the optimization
% % Initial State :
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
q_AIRcathinp  = 17.2731441172150; % [L/min]     methane flow rate
Iinp          = 19;               % [A]         current

u_0 = [q_CH4inp,q_AIRcathinp,Iinp];

% Initial guess for the optimization problem
u0 = [u_0 T_0];

%% Constraints 
Ps_el = 100;

% ub = [27.25,272.53/0.21, 30, 600*ones(1,1)+273.15 800*ones(1,4)+273.15 1600+273.15 1474 1474];
% lb = [1.36E-03,0.01/0.21,0, 450*ones(1,6)+273.15   200+273.15 200+273.15];

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


%% Optimization 
myProblem.TC.TimeConstantOff = 1;
tic
[u_f,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,[],T_in,A,B,lb,ub,SOFC_data_nominal,myProblem.OPT.options);
[eff,myc,myceq] = computeall(u_f,[],T_in,SOFC_data_nominal);
toc

%%
[tempsteady_plant] = OPTIM_SteadyState(u_f(1:3),u_f(4:end),T_in,SOFC_data_nominal);
[dx_global,Ucell,Pel,FU,L_air,eta_el,eta_sys] = fPrimeMyProject(0,tempsteady_plant,u_f(1:3),T_in,SOFC_data_nominal);
%%
nCH4in  = u_f(1)
nAIRin  = u_f(2)
I       = u_f(3)
FU      = kc*u_f(3)/u_f(1)
Lair    = u_f(2)/2/u_f(1)
T_reformer    = u_f(4) - 273.15

T_SOFC = mean(u_f(5:8)) - 273.15

Tburn      = u_f(9) - 273.15
T_HEXfuel  = u_f(10) - 273.15
T_HEXair   = u_f(11) - 273.15

