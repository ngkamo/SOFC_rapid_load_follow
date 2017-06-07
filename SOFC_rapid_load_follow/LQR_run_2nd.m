clear all;
close all;
clc;
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
initial = 1;
target = 2;
Pel_opt = [90 100];
Ps_el = [Pel_opt(1)*ones(1,5) Pel_opt(2)*ones(1,5)]; % Power profile

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
for i = 1:2
    Ps_el = Pel_opt(i);

    [u_f,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,[],T_in,A,B,lb,ub,SOFC_data_nominal,myProblem.OPT.options);
    [dx,uc_opt,~,~,~,~,eff_opt(i)] = fPrimeMyProject(0,u_f(4:end),u_f(1:3),T_in,SOFC_data_nominal);
    u_opt = [u_opt; u_f];
    dx_opt = [dx_opt, dx];
    Ucell_opt = [Ucell_opt, uc_opt];
end

%% Testing delta stability
% delta = logspace(-8,0,50);
% index = 3;
% differential = [];
% for i = 1:size(delta,2)
%     u_grad = u_f;
%     u_grad(index) = u_grad(index) + delta(i);
% 
%     [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,eff_test] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
% 
%     differential(1,i) = (eff_test-eff_opt(target))/delta(i);
% end

% loglog(delta,abs(differential));

%% ------------------------------
% Linearizing the system
% -------------------------------
[A,B,C,D] = MPC_linearization_mat_LQR(u_f,dx_opt(:,target),Pel_opt(target),Ucell_opt(target),eff_opt(target),T_in,SOFC_data_nominal);

ns = size(A,1);  % number of states
ni = size(B,2);  % number of inputs
no = size(C,1);  % number of outputs

%% Continuous time state-space model
sys = ss(A,B,C,D);

Q = C'*C + 10*eye(8);
R = 3*eye(3);
K = lqr(sys,Q,R);

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;

%% Discrete time state-space model
Ts = 10;
sys_d = c2d(sys,Ts,'zoh');
Q = sys_d.C'*sys_d.C + 0.1*eye(8);
R = 100*eye(3);
Kd = dlqr(sys_d.A,sys_d.B,Q,R);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

%% Initial/target states
H = [Cd Dd; Ad-eye(ns) Bd];
ytarget = [0; 0; zeros(8,1)];
var_target = H\ytarget;

yinitial = [Pel_opt(initial)-Pel_opt(target); 0; zeros(ns+ni-no-1,1)];
var_initial = H\yinitial;

%% Simulation with linearized model
x_ss = u_opt(target,4:end)';  % states corresponding to steady state
u_ss = u_opt(target,1:3)';
x_initial = u_opt(initial,4:end)';
x_target = x_ss;

Nsim = 200;

% initial states
x_hat_initial = (x_initial-x_ss);
x_hat_target = (x_target-x_ss);
x_hat = var_initial(1:8);
var_tracking = var_initial;
x_lin = var_initial(1:8)+x_ss;
u_hat = [];
y_hat = [];

for i = 1:Nsim
    if i > 75
        var_tracking = var_target;
    else
        var_tracking = var_initial;
    end
    u_hat(:,i) = -Kd*(x_hat(:,i)-var_tracking(1:8)) + var_tracking(9:11);
    u_min = 0.12;
    
    if u_hat(1,i)+u_opt(target,1)< u_min
        u_hat(1,i) = -u_opt(target,1) + u_min;
    end

    x_hat(:,i+1) = Ad*x_hat(:,i) + Bd*u_hat(:,i);
    x_lin(:,i+1) = x_hat(:,i+1) + x_ss;
    y_hat(:,i)   = Cd*x_hat(:,i) + Dd*u_hat(:,i);
end
y_hat = [y_hat(:,1) y_hat];

%% Simulation with non-linearized model
myProblem.TC.TimeConstantOff = 0;
u_ss = u_opt(target,1:3)';
x_hat_non = var_initial(1:8);
x_nonlin = x_hat_non + x_ss;

opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

% sol = ode15s(@(t,x) SIMU_ode_func(t,x,T_in,SOFC_data_nominal,Kd,x_ss,u_ss,var_target) , time, x0, opts);
for i = 1:Nsim
    if i > 75
        var_tracking = var_target;
    else
        var_tracking = var_initial;
    end
    u_nonlin(:,i) = -Kd*(x_hat_non(:,i)-var_tracking(1:8)) + var_tracking(9:11) + u_ss;
    u_min = 0.012;
    if u_nonlin(1,i) < u_min
        u_nonlin(1,i) = u_min ;
    end
    
    sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_nonlin(:,i),T_in,SOFC_data_nominal), [0 Ts], x_nonlin(:,i), opts);
    x_nonlin(:,i+1) = sol.y(:,end);
    [~,U_nonlin(i),P_nonlin(i),FU_nonlin(i),L_nonlin(i),~,eta_nonlin(i)] = fPrimeMyProject(0,x_nonlin(:,i+1),u_nonlin(:,i),T_in,SOFC_data_nominal);
    x_hat_non(:,i+1) = x_nonlin(:,i+1) - x_ss;
end

%% 
close all
time = Ts/60*(0:1:Nsim); % [s]
set(0,'DefaultFigureWindowStyle','docked')

figure('Name','temperature')
[ts,xs] = stairs(time,x_lin');
plot(ts,xs(:,6)',time,x_nonlin(6,:));
ylabel('Temperatures [K]')
xlabel('Time [min]');
legend('linear','nonlinear')
grid on

figure('Name','power')
[ts,ys] = stairs(time(1:end),y_hat(1,:)+Pel_opt(target)');
plot(ts,ys,time(2:end),P_nonlin)
ylabel('Power [W]')
xlabel('Time [min]')
title('Power')
legend('linear','nonlinear')
grid on
print('power','-dpng')

figure('Name','voltage')
[ts,ys] = stairs(time(1:end),-y_hat(2,:)+Ucell_opt(target));
plot(ts,ys, time(2:end),U_nonlin)
ylabel('Voltage [V]')
xlabel('Time [min]')
title('Voltage')
legend('linear','nonlinear')
grid on
print('voltage','-dpng')
 
% figure('Name','efficiency')
% [ts,ys] = stairs(time(1:end),y_hat(3,:)+eff_opt(target));
% plot(ts,ys, time(2:end),eta_nonlin)
% ylabel('Efficiency \eta[-]')
% xlabel('Time [min]')
% title('Efficiency')
% legend('linear','nonlinear')
% grid on
% print('efficiency','-dpng')

figure('Name','CH4')
[ts,us] = stairs(time(1:end-1),u_hat'+u_f(1:3));
plot(ts(:,1),us(:,1),time(2:end),u_nonlin(1,:))
ylabel('Flow rate CH_4 [L/min]')
xlabel('Time [min]')
title('Flow rate CH_4')
legend('linear','nonlinear')
grid on
print('ch4','-dpng')

figure('Name','air')
plot(ts(:,1),us(:,2),time(2:end),u_nonlin(2,:))
ylabel('Flow rate air [L/min]')
xlabel('Time [min]')
title('Flow rate air')
legend('linear','nonlinear')
grid on
print('air','-dpng')

figure('Name','current')
plot(ts(:,1),us(:,3),time(2:end),u_nonlin(3,:))
ylabel('Current [A]')
xlabel('Time [min]')
title('Current')
legend('linear','nonlinear')
grid on
print('current','-dpng')