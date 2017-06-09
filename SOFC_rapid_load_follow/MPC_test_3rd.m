clear all
close all
yalmip('clear')
addpath('REFORMER','SOFC','COR','HEX','BURN');

%%
global myProblem Ps_el it u_previous

% Number of cells :
N_c = 6;

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);
SOFC_data_plant = data_SOFC_plant(prjname_SOFC,N_c);
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

%% ------------------------------
% Linearizing the system
% -------------------------------
[Ac,Bc,Cc,Dc] = MPC_linearization_mat(u_f,dx_opt(:,target),Pel_opt(target),Ucell_opt(target),T_in,SOFC_data_nominal);

%% Continuous time state-space model
sys = ss(Ac,Bc,Cc,Dc);

Q = Cc'*Cc+0.01*eye(8);
R = 3*eye(3);
K = lqr(sys,Q,R);

%% Discrete time state-space model
Ts = 10;
sys_d = c2d(sys,Ts,'zoh');
% Q = sys_d.C'*sys_d.C + 0.1*eye(8);
% R = 10*eye(3);
% Kd = dlqr(sys_d.A,sys_d.B,Q,R);
A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

%% Steady-state target computation
ns = 8;
ni = 3;
no = 2;

% Variables
x_target = sdpvar(ns,1,'full');
u_target = sdpvar(ni,1,'full');
rs  = sdpvar(no,1,'full');

% r_target = [10; 0; -0.01];
r_target = [[0; 0].*ones(2,50), zeros(2,300)]; % profile demand

% Constraints and objective
% ub = [27.25; 272.53/0.21; 30];
% ub = [27.25; 120; 30];
% lb = [1.36E-03; 0.01/0.21; 0];
ub = [27.25; 120; 30];
lb = [0; 0; 0];
Lair_upper = 10; 
Lair_lower = 3;
FU_upper   = 0.7;

Rcon = 8.314462175;    % [J/K/mol] 
kc = (6e+4)*N_c*Rcon*T_st/(8*P*SOFC_data_nominal.cst.F);
M   = [2*Lair_lower,  -0.21,  0;
        -2*Lair_upper,   0.21,  0;
        -FU_upper,          0,  kc];
m   = [0;0;0];

x_ss = u_opt(target,4:end)';
u_ss = u_opt(target,1:3)';
y_ss = [100; 0.7];

H = [C(1:2,:) D(1:2,:); A-eye(ns) B];
% H\[10; 0; zeros(8,1)]

con = [];
obj = 0;
R = ones(3);

con = con + ( x_target == A*x_target + B*u_target );
con = con + ( rs == C(1:2,:)*x_target + D(1:2,:)*u_target );
con = con + ( M*(u_target+u_ss) <= m );
% con = con + ( lb <= (u_target)  <= ub);
% con = con + (0.5 <= u_target(2)<= 100);
ob = obj + x_target'*Q*x_target + u_target'*R*u_target;

% Parameters
parameters_in = {rs};
parameters_out = {x_target, u_target};

% Compile the matrices
steady_target = optimizer(con, obj, [], parameters_in, parameters_out);
[solutions_initial,infeasible] = steady_target{[-10; 0]};
[solutions_target,infeasible] = steady_target{[0; 0]};
solutions_initial{1}
solutions_initial{2}
solutions_initial{1}'*Q*solutions_initial{1} + solutions_initial{2}'*R*solutions_initial{2}
C(1:2,:)*solutions_initial{1} + D(1:2,:)*solutions_initial{2}

solutions_initial = H\[Pel_opt(initial)-Pel_opt(target); 0; zeros(ns+ni-no-1,1)];
solutions_target = H\[0; 0; zeros(ns+ni-no-1,1)];

%% MPC tracking
Q = 1*eye(ns);
R = 100*eye(ni);

N = 40;  % horizon length

% Variables
x_hat  = sdpvar(ns,N,'full');
x_target = sdpvar(ns,1,'full');
u_hat  = sdpvar(ni,N,'full');
u_target = sdpvar(ni,1,'full');

% Define constrainst and objective
con = [];
obj = 0;

for i = 1:N-1
    con = con + ( x_hat(:,i+1) == A*x_hat(:,i) + B*u_hat(:,i) );
    con = con + ( M*(u_hat(:,i)+u_ss) <= m );
    con = con + ( lb <= u_hat(:,i)+u_ss <= ub );
    
    y(:,i) = C*x_hat(:,i) + D*u_hat(:,i) + y_ss;
    con = con + ( y(2,i) >= 0.7 );
    
    obj = obj + (x_hat(:,i)-x_target)'*Q*(x_hat(:,i)-x_target)...
              + (u_hat(:,i)-u_target)'*R*(u_hat(:,i)-u_target);
end

slew_rate = repmat([0.0005 1 2],N-1,1);
test = diff(u_hat,1,2);
con = con + ( -slew_rate <= diff(u_hat,1,2)' <= slew_rate);
con = con + ( M*(u_hat(:,N)-u_target+u_ss) <= m );
con = con + ( lb <= u_hat(:,N)-u_target+u_ss <= ub );

obj = obj + (x_hat(:,i)-x_target)'*Q*(x_hat(:,i)-x_target)...
          + (u_hat(:,i)-u_target)'*R*(u_hat(:,i)-u_target);
      



% Parameters
parameters_in = {x_hat(:,1), x_target, u_target};
solutions_out = {u_hat};

% Compile the matrices
controller = optimizer(con, obj, [], parameters_in, solutions_out);

%% Simulation linear model
myProblem.TC.TimeConstantOff = 0;
Nsim = size(r_target,2);  % length of the simulation
Nsim = 180;
x_hat = [];
y = [];
u_lin = [];
x_hat(:,1) = solutions_initial(1:8);%(u_opt(2,4:end)'-x_ss);

for j = 1:Nsim
    if j < 60
        x_target = solutions_initial(1:8);
        u_target = solutions_initial(9:11);
    else
        x_target = solutions_target(1:8);
        u_target = solutions_target(9:11);
    end
    [solutions, infeasible] = controller{{x_hat(:,j), x_target, u_target}};
    u_lin(:,j) = solutions(:,1);
    
    x_hat(:,j+1) = A*x_hat(:,j) + B*u_lin(:,j);
    y(:,j)   = C(1:2,:)*x_hat(:,j) + D(1:2,:)*u_lin(:,j) + y_ss;
end

%% Simulation non-linear model
myProblem.TC.TimeConstantOff = 0;
x_ss = u_opt(target,4:end)';
u_ss = u_opt(target,1:3)';
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

x_initial = u_opt(initial,4:end)'-u_opt(target,4:end)';
u_initial = u_opt(initial,1:3)'-u_opt(initial,1:3)';

x_hat_nonlin = solutions_initial(1:8);
% x_hat(:,1) = (u_opt(2,4:end)'-x_ss);
x_nonlin = x_hat_nonlin + x_ss;

for j = 1:Nsim
    if j < 60
        x_setpoint = solutions_initial(1:8);
        u_target = solutions_initial(9:11);
    else
        x_setpoint = solutions_target(1:8);
        u_target = solutions_target(9:11);
    end
    [u_hat, infeasible] = controller{{x_hat_nonlin(:,j), x_setpoint, u_target}};
    u_nonlin(:,j) = u_hat(:,1)+u_ss;
    
    sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_nonlin(:,j),T_in,SOFC_data_plant), [0 Ts], x_nonlin(:,j), opts);
    x_nonlin(:,j+1) = sol.y(:,end);
    x_hat_nonlin(:,j+1) = x_nonlin(:,j+1) - x_ss;
    [~,U_nonlin(j),P_nonlin(j),FU_nonlin(j),L_nonlin(j),~,eta_nonlin(j)] = fPrimeMyProject(0,x_nonlin(:,j+1),u_nonlin(:,j),T_in,SOFC_data_plant);
end

%%
close all
timestep = Ts*[0:Nsim]/60;
set(0,'DefaultFigureWindowStyle','docked')
figure('Name','power')
[ts,ys] = stairs(timestep(1:end-1),y(1,:)');
plot(ts,ys,timestep(1:end-1),P_nonlin)
xlabel('Time [min]')
ylabel('Power [W]')
title('Power [W]')
legend('linear','nonlinear')
grid on

figure('Name','voltage')
[ts,ys] = stairs(timestep(1:end-1),y(2,:));
plot(ts,ys, timestep(1:end-1),U_nonlin)
xlabel('Time [min]')
ylabel('Voltage [V]')
title('Voltage [V]')
legend('linear','non linear')
grid on

% figure('Name','efficiency')
% [ts,ys] = stairs(timestep(1:end-1),y(3,:));
% plot(ts,ys, timestep(1:end-1),eta_nonlin);
% xlabel('Time [min]')
% ylabel('Efficiency [-]')
% title('Efficiency')
% legend('linear','nonlinear')
% grid on

figure('Name','methane')
[ts,ys] = stairs(timestep(1:end-1),u_lin(1,:)+u_opt(target,1));
plot(ts,ys, timestep(1:end-1),u_nonlin(1,:))
xlabel('Time [min]')
ylabel('q_{CH4} [L/min]')
title('Input methane [L/min]')
legend('linear','non linear')
grid on

figure('Name','air')
[ts,ys] = stairs(timestep(1:end-1),u_lin(2,:)+u_opt(target,2));
plot(ts,ys, timestep(1:end-1),u_nonlin(2,:))
xlabel('Time [min]')
ylabel('q_{air} [L/min]')
title('Input air flow [L/min]')
legend('linear','non linear')
grid on

figure('Name','current')
[ts,ys] = stairs(timestep(1:end-1),u_lin(3,:)+u_opt(target,3));
plot(ts,ys, timestep(1:end-1),u_nonlin(3,:));
xlabel('Time [min]')
ylabel('current I [A]')
title('Current [A]')
legend('linear','non linear')
grid on

figure('Name','temperatures')
[ts,ys] = stairs(timestep(1:end),x_hat(6,:));
plot(ts,ys,timestep,x_hat_nonlin(6,:));
xlabel('Time [min]')
ylabel('Temperature [K]')
title('Temperature [K]')
legend('linear','non linear')
grid on