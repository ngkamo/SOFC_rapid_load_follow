function [controller, target_config,U_nom,P_nom,eff_nom] = MPC_controller(optimal_var,Ts, x_ss,u_ss,y_ss,Pel_target,T_in,SOFC_data_nominal,M,m,lb,ub,slew_rate)
ns = 8; % # states
ni = 3; % # inputs
no = 2; % # outputs

%% Linearizing the system

[dotx,U_nom,P_nom,~,~,~,eff_nom] = fPrimeMyProject(0,x_ss,u_ss,T_in,SOFC_data_nominal);
[Ac,Bc,Cc,Dc] = MPC_linearization_mat(optimal_var,dotx,Pel_target,U_nom,T_in,SOFC_data_nominal);

% Discrete time state-space model
sys = ss(Ac,Bc,Cc,Dc);

sys_d = c2d(sys,Ts,'zoh');

A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

% Linearized steady states and input @ power setpoint
H = [C D; A-eye(ns) B];

target_config = H\[0; 0; zeros(ns+ni-no-1,1)];

%% MPC controller
Q = 1*eye(ns);
% R_scaling = [ 1000/0.3   1/15   1/20 ];
% R = 100*diag(R_scaling);
R = 100*eye(ni);

N = 40;  % horizon length

% Variables
x_hat  = sdpvar(ns,N,'full');
x_target = sdpvar(ns,1,'full');
u_hat  = sdpvar(ni,N,'full');
% u_nonlin = sdpvar(ni,N,'full');
u_target = sdpvar(ni,1,'full');
u_prev = sdpvar(ni,1,'full');
y = sdpvar(no,N,'full');

% Define constrainst and objective
con = [];
obj = 0;

con = con + ( -slew_rate <= u_hat(:,1)+u_ss-u_prev <= slew_rate );

for i = 1:N-1
    con = con + ( x_hat(:,i+1) == A*x_hat(:,i) + B*u_hat(:,i) );
    con = con + ( M*(u_hat(:,i)+u_ss) <= m );
    con = con + ( lb <= u_hat(:,i)+u_ss <= ub );
    
    con = con + ( -slew_rate <= u_hat(:,i+1)-u_hat(:,i) <= slew_rate );
    
    con = con + ( y(:,i) == C*x_hat(:,i) + D*u_hat(:,i) + y_ss );
    con = con + ( y(2,i) >= 0.7 );
    
    obj = obj + (x_hat(:,i)-x_target)'*Q*(x_hat(:,i)-x_target)...
              + (u_hat(:,i)-u_target)'*R*(u_hat(:,i)-u_target);
end

con = con + ( M*(u_hat(:,N)-u_target+u_ss) <= m );
con = con + ( lb <= u_hat(:,N)-u_target+u_ss <= ub );
obj = obj + (x_hat(:,i)-x_target)'*Q*(x_hat(:,i)-x_target)...
          + (u_hat(:,i)-u_target)'*R*(u_hat(:,i)-u_target);

% Parameters
parameters_in = {x_hat(:,1), x_target, u_target, u_prev};
solutions_out = {u_hat};

% Compile the matrices
controller = optimizer(con, obj, [], parameters_in, solutions_out);

end