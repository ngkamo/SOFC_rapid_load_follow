function [states_plant, input_plant, last_config, P_plant, U_plant, FU_plant, Lair_plant, eff_plant, output_nompla, modif] = MPC_epsilon(optimal_var,modif,Pel_target,Ucell_target,T_in,SOFC_data_nominal,SOFC_data_plant,M,m,lb,ub,initial_config)
% global myProblem
yalmip('clear')

Ts = 10;  % Time sample

x_ss = optimal_var(4:end);
u_ss = optimal_var(1:3);
y_ss = [Pel_target; Ucell_target];

% Reformating constraints
M(:,4:end) = [];
lb = lb';
ub = ub';
lb(4:end)  = [];
ub(4:end)  = [];

slew_rate = [0.001; 0.05; 0.3];

%% MPC controller
[controller,target_config,U_nom,P_nom,eff_nom] = MPC_controller(optimal_var,Ts,x_ss,u_ss,y_ss,Pel_target,T_in,SOFC_data_nominal,M,m,lb,ub,slew_rate);

%% Simulation non-linear plant model
Nsim = 300;
tol  = 1e-4;
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

x_setpoint = target_config(1:8);
u_target = target_config(9:11);
x_plant = initial_config(4:end);
xhat_plant = x_plant-x_ss;

u_previous = initial_config(1:3);

for j = 1:Nsim
    [u_hat, infeasible] = controller{{xhat_plant(:,j), x_setpoint, u_target, u_previous}};
    u_plant(:,j) = u_hat(:,1)+u_ss;
    u_plant(:,j) = MPC_apply_slewrate(u_plant(:,j),u_previous,slew_rate);
    u_previous = u_plant(:,j);
    
    sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_plant(:,j),T_in,SOFC_data_plant), [0 Ts], x_plant(:,j), opts);
    x_plant(:,j+1) = sol.y(:,end);
    xhat_plant(:,j+1) = x_plant(:,j+1) - x_ss;
    [dx_states,U_plant(j),P_plant(j),FU_plant(j),Lair_plant(j),~,eff_plant(j)] = fPrimeMyProject(0,x_plant(:,j+1),u_plant(:,j),T_in,SOFC_data_plant);
    
    check_steady_plant = norm(dx_states);
    if check_steady_plant < tol
        break;
    end
end

input_plant  = u_plant;
states_plant = x_plant(:,1:end-1);
last_config  = [u_plant(:,end); x_plant(:,end-1)];

%% Computing the constraint adaptation modifiers
Kca = [0 0.9 0.9];
modif(2,4) = (1-Kca(2))*modif(2,4) + Kca(2)*(U_plant(end)-U_nom(end));    % Cell voltage
modif(3,4) = (1-Kca(3))*modif(3,4) + Kca(3)*(P_plant(end)-P_nom(end));        % Power set point

output_nompla = [P_nom, U_nom, eff_nom, P_plant(end), U_plant(end), eff_plant(end)];
end