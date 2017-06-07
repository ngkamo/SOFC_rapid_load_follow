clear all;
close all;
clc;
warning on
addpath('REFORMER','SOFC','COR','HEX','BURN');

%%
global myProblem Ps_el it u_previous

% Number of cells :
N_c = 6;

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c)
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
q_AIRcathinp  = 17.2531441172150; % [L/min]     methane flow rate
Iinp          = 19;               % [A]         current

u_0 = [q_CH4inp,q_AIRcathinp,Iinp];

% Initial guess for the optimization problem
u0 = [u_0 T_0];

%% Constraints 
target = 3;
Pel_opt = [80 100 90];
Ps_el = [Pel_opt(1)*ones(1,5) Pel_opt(2)*ones(1,5) Pel_opt(3)*ones(1,5)]; % Power profile

% ub = [27.25,272.53/0.21, 30, 600*ones(1,1)+273.15 800*ones(1,4)+273.15 1600+273.15 1474 1474];
% lb = [1.36E-03,0.01/0.21, 0, 450*ones(1,6)+273.15   200+273.15 200+273.15];
ub = [27.25,272.53/0.21, 30, Inf*ones(1,8)];
lb = [1.36E-03,0.01/0.21, 0, -Inf*ones(1,8)];

Aeq = [];
beq = [];

Lair_upper = 10; 
Lair_lower = 3;
FU_upper   = 0.7;

kc = (6e+4)*N_c*R*T_st/(8*P*SOFC_data_nominal.cst.F);

A   = [2*Lair_lower,  -0.21,  0,  zeros(1,8);
      -2*Lair_upper,   0.21,  0,  zeros(1,8);
      -FU_upper,          0,  kc, zeros(1,8)];
B   =  [0;0;0];

%% -----------------------------------------------
% Optimal values of the plant at power setpoints
% ------------------------------------------------
myProblem.TC.TimeConstantOff = 0;
u_opt = [];
eff_opt = [];
dx_opt = [];
Ucell_opt = [];

% operating power
for i = 1:3
    Ps_el = Pel_opt(i);

    [u_f,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,[],T_in,A,B,lb,ub,SOFC_data_nominal,myProblem.OPT.options);
    [dx,uc_opt,~,~,~,~,eff_opt(i)] = fPrimeMyProject(0,u_f(4:end),u_f(1:3),T_in,SOFC_data_nominal);
    u_opt = [u_opt; u_f];
    dx_opt = [dx_opt, dx];
    Ucell_opt = [Ucell_opt, uc_opt];
end

%% Testing delta stability
delta = logspace(-8,0,50);
index = 3;
differential = [];
for i = 1:size(delta,2)
    u_grad = u_f;
    u_grad(index) = u_grad(index) + delta(i);

    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,eff_test] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);

    differential(1,i) = (eff_test-eff_opt(target))/delta(i);
end

% loglog(delta,abs(differential));

%% ------------------------------
% Linearizing the system
% -------------------------------
% A & B
u_f = u_opt(target,:);
deltaA = 0.001*ones(1,8);
differentialA = [];
component = 8;
A = zeros(8,8);
B = zeros(8,3);

for i = 1:length(deltaA)
    u_grad = u_f;
    u_grad(3+i) = u_grad(3+i) + deltaA(i);
    
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    differentialA(:,1) = (dx_test-dx_opt(:,target))/deltaA(i);
    
    A(:,i) = differentialA;
end

deltaB = [1e-6 1e-5 1e-5];
differentialB = [];
for i = 1:length(deltaB)
    u_grad = u_f;
    u_grad(i) = u_grad(i) + deltaB(i);
    
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    differentialB(:,1) = (dx_test-dx_opt(:,target))/deltaB(i);
    
    B(:,i) = differentialB;
end

% C & D
deltaCPel = 1e-2*ones(1,8);
deltaCUc  = 1e-2*ones(1,8);
deltaCEff = 1e-4*ones(1,8);
differentialCPel = [];
differentialCUc  = [];
differentialCEff = [];
C = zeros(3,8);

for i = 1:length(deltaCPel)
    u_grad = u_f;
    u_grad(3+i) = u_grad(3+i) + deltaCPel(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    C(1,i) = (Pel_test-Pel_opt(target))/deltaCPel(i);
    C(2,i) = (Ucell_test-Ucell_opt(target))/deltaCUc(i);
    
    u_grad = u_f;
    u_grad(3+i) = u_grad(3+i) + deltaCPel(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,eff_test] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    C(3,i) = (eff_test-eff_opt(target))/deltaCEff(i);
end

deltaDPel = [1e-8 1e-6 1e-6];
deltaDUc  = [1e-8 1e-2 1e-3];
deltaDEff = [1e-7 1e-2 1e-2];
differentialD = [];
D = zeros(3,3);
for i = 1:length(deltaDPel)
    u_grad = u_f;
    u_grad(i) = u_grad(i) + deltaDPel(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    D(1,i) = (Pel_test-Pel_opt(target))/deltaDPel(i);
    
    u_grad = u_f;
    u_grad(i) = u_grad(i) + deltaDUc(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    D(2,i) = (Ucell_test-Ucell_opt(target))/deltaDUc(i);
    
    u_grad = u_f;
    u_grad(i) = u_grad(i) + deltaDUc(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,eff_test] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
    D(3,i) = (eff_test-eff_opt(target))/deltaDEff(i);
end

%% Continuous time state-space model
sys = ss(A,B,C,D);

Q = C'*C+0.01*eye(8);
R = 3*eye(3);
K = lqr(sys,Q,R);

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;

%% Discrete time state-space model
Ts = 1;
sys_d = c2d(sys,Ts,'zoh');
Q = sys_d.C'*sys_d.C + 1*eye(8);
R = 3*eye(3);
Kd = dlqr(sys_d.A,sys_d.B,Q,R);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

H = [Cd , Dd; Ad-eye(8) Bd];
% ytarget = [-10; 0; -0.0029; zeros(8,1)];
ytarget = [0; 0; 0; zeros(8,1)];
var_target = H\ytarget;

%% Simulation
Nsim = 200;
% ytest = 0 + C*(u_opt(1,4:end)'-u_f(4:end)') + D*(u_opt(2,1:3)'-u_f(1:3)');
% x = -[10.8179  -14.1690  -14.0954  -14.0930  -14.0878  -30.8718  -18.8475  -18.0659]';
% x = -[var_target(1:8)];
x = [10 0 0 0 0 0 0 0]';
% x = zeros(8,1);
y0 = Cd*x(:,1)-Dd*[var_target(9:11)];
u = [];
y = [];

for i = 1:Nsim
    u(:,i) = -Kd*(x(:,i)-var_target(1:8)) + var_target(9:11);
%     u(:,i) = -Kd*x(:,i);
    x(:,i+1) = Ad*x(:,i) + Bd*u(:,i);
    y(:,i) = Cd*x(:,i) + Dd*u(:,i);
end
y = [y(:,1) y];

%%
t = Ts*(0:1:Nsim);
figure(1)
[ts,ys] = stairs(t(1:end),y(1,:)+Pel_opt(target));
plot(ts,ys)
ylabel('Power [W]')
xlabel('Time [s]')
title('Power')
grid on
print('power','-dpng')

figure(2)
[ts,ys] = stairs(t(1:end),-y(2,:)+Ucell_opt(target));
plot(ts,ys)
ylabel('Voltage [V]')
xlabel('Time [s]')
title('Voltage')
grid on
print('voltage','-dpng')

figure(3)
[ts,ys] = stairs(t(1:end),y(3,:)+eff_opt(target));
plot(ts,ys)
ylabel('Efficiency \eta[-]')
xlabel('Time [s]')
title('Efficiency')
grid on
print('efficiency','-dpng')

figure(4)
[ts,xs] = stairs(t,x');
plot(ts(:,1),xs(:,2))
ylabel('Temperatures [K]')
xlabel('Time [s]');
grid on

figure(5)
[ts,us] = stairs(t(2:end),u'+u_f(1:3));
plot(ts(:,1),us(:,1))
ylabel('Flow rate CH_4 [L/min]')
xlabel('Time [s]')
title('Flow rate CH_4')
grid on
print('ch4','-dpng')

figure(6)
plot(ts(:,1),us(:,2))
ylabel('Flow rate air [L/min]')
xlabel('Time [s]')
title('Flow rate air')
grid on
print('air','-dpng')

figure(7)
plot(ts(:,1),us(:,3))
ylabel('Current [A]')
xlabel('Time [s]')
title('Current')
grid on
print('current','-dpng')