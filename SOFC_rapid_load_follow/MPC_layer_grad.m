% MPC layer
function [output_plant, input_plant, last_config, P_plant, U_plant, eff_plant, modif] = MPC_layer(optimal_var,modif,Pel_target,Ucell_target,T_in,SOFC_data_nominal,SOFC_data_plant,M,m,lb,ub,initial_config)
global myProblem

yalmip('clear')

ns = 8; % # states
ni = 3; % # inputs
no = 2; % # outputs

x_ss = optimal_var(ni+1:end);
u_ss = optimal_var(1:3);
y_ss = [Pel_target; Ucell_target];

% Reformating constraints
M(:,4:end) = [];
lb = lb';
ub = ub';
lb(4:end)  = [];
ub(4:end)  = [];

%% ------------------------------
% Linearizing the system
% -------------------------------

[dotx,Ucell_opt,~,~,~,~,~] = fPrimeMyProject(0,optimal_var(4:end),optimal_var(1:3),T_in,SOFC_data_nominal);
[Ac,Bc,Cc,Dc] = MPC_linearization_mat(optimal_var,dotx,Pel_target,Ucell_opt,T_in,SOFC_data_nominal);

% Discrete time state-space model
sys = ss(Ac,Bc,Cc,Dc);

Ts = 10;  % Time sample
sys_d = c2d(sys,Ts,'zoh');

A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

% Linearized steady states and input @ power setpoint
H = [C D; A-eye(ns) B];

target_config = H\[0; 0; zeros(ns+ni-no-1,1)];

%% MPC tracking
Q = 1*eye(ns);
R = 100*eye(ni);

N = 20;  % horizon length

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

con = con + ( M*(u_hat(:,N)-u_target+u_ss) <= m );
con = con + ( lb <= u_hat(:,N)-u_target+u_ss <= ub );
obj = obj + (x_hat(:,i)-x_target)'*Q*(x_hat(:,i)-x_target)...
          + (u_hat(:,i)-u_target)'*R*(u_hat(:,i)-u_target);

% Parameters
parameters_in = {x_hat(:,1), x_target, u_target};
solutions_out = {u_hat};

% Compile the matrices
controller = optimizer(con, obj, [], parameters_in, solutions_out);

%% Simulation non-linear plant model
Nsim = 200;
tol  = 9e-3;
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

x_setpoint = target_config(1:8);
u_target = target_config(9:11);

xhat_nonlin = initial_config(4:end);
x_nonlin = initial_config(4:end);
xhat_nonlin = x_nonlin-x_ss;

for j = 1:Nsim
    [u_hat, infeasible] = controller{{xhat_nonlin(:,j), x_setpoint, u_target}};
    u_nonlin(:,j) = u_hat(:,1)+u_ss;
    
    sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_nonlin(:,j),T_in,SOFC_data_plant), [0 Ts], x_nonlin(:,j), opts);
    x_nonlin(:,j+1) = sol.y(:,end);
    xhat_nonlin(:,j+1) = x_nonlin(:,j+1) - x_ss;
    [~,U_plant(j),P_plant(j),~,~,~,eff_plant(j)] = fPrimeMyProject(0,x_nonlin(:,j+1),u_nonlin(:,j),T_in,SOFC_data_plant);
    
%     check_ss = norm(dx_states);
    check_steady_plant = norm(diff(x_nonlin(:,j:j+1),1,2));
    if check_steady_plant < tol
        break;
    end
end

input_plant  = u_nonlin;
output_plant = x_nonlin(:,1:end-1);
last_config = [u_nonlin(:,end); x_nonlin(:,end-1)];

for i = 1:3
    for j = 1:Nsim
        [u_hat, infeasible] = controller{{xhat_nonlin(:,j), x_setpoint, u_target}};
        u_nonlin(:,j) = u_hat(:,1)+u_ss;

        sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_nonlin(:,j),T_in,SOFC_data_plant), [0 Ts], x_nonlin(:,j), opts);
        x_nonlin(:,j+1) = sol.y(:,end);
        xhat_nonlin(:,j+1) = x_nonlin(:,j+1) - x_ss;
        [~,U_plant(j),P_plant(j),~,~,~,eff_plant(j)] = fPrimeMyProject(0,x_nonlin(:,j+1),u_nonlin(:,j),T_in,SOFC_data_plant);
        
        check_steady_plant = norm(diff(x_nonlin(:,j:j+1),1,2));
        if check_steady_plant < tol
            break;
        end
    end
    
end


%% Simulation non-linear nominal model
Nsim = 200;
% opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

x_setpoint = target_config(1:8);
u_target = target_config(9:11);

xhat_nom = initial_config(4:end);
x_nom = initial_config(4:end);
xhat_nom = x_nom-x_ss;

for j = 1:Nsim
    [u_hat, infeasible] = controller{{xhat_nom(:,j), x_setpoint, u_target}};
    u_nonlin(:,j) = u_hat(:,1)+u_ss;
    
    sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_nonlin(:,j),T_in,SOFC_data_nominal), [0 Ts], x_nom(:,j), opts);
    x_nom(:,j+1) = sol.y(:,end);
    xhat_nom(:,j+1) = x_nom(:,j+1) - x_ss;
    [~,U_nom(j),P_nom(j),~,~,~,~] = fPrimeMyProject(0,x_nom(:,j+1),u_nonlin(:,j),T_in,SOFC_data_nominal);
    
%     check_ss = norm(dx_states);
    check_steady_nom = norm(diff(x_nom(:,j:j+1),1,2));
    if check_steady_nom < tol
        break;
    end
end

% [tempsteady_nominal] = OPTIM_SteadyState(input_plant(:,end)',output_plant(:,end)',T_in,SOFC_data_nominal);
% [~,U_nom,P_nom,~,~,~,~] = fPrimeMyProject(0,tempsteady_nominal,u_nonlin(:,j),T_in,SOFC_data_nominal);


%% Computing the modifiers
Kca = [0 1 1];
modif(2,4) = (1-Kca(2))*modif(2,4) + Kca(2)*(U_plant(end)-U_nom(end));    % Cell voltage
modif(3,4) = (1-Kca(3))*modif(3,4) + Kca(3)*(P_plant(end)-P_nom(end));        % Power set point

end