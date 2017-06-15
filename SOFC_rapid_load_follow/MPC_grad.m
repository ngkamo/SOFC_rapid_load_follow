function [output_plant, input_plant, last_config, P_plant_hist, U_plant_hist, eff_plant_hist, modif] = MPC_grad(optimal_var,output_nompla,modif,Pel_target,Ucell_target,T_0,T_in,SOFC_data_nominal,SOFC_data_plant,M,m,lb,ub,initial_config,last_plant_output)
Ts = 10;  % Time sample

x_ss = optimal_var(4:end);
u_ss = optimal_var(1:3);
input_plant  =  [];
output_plant = [];
P_plant_hist = [];
eff_plant_hist = [];
U_plant_hist = [];
% y_ss = [Pel_target; Ucell_target];

P_nom1 = output_nompla(1);
U_nom1 = output_nompla(2);
eff_nom1 = output_nompla(3);
P_pla1 = output_nompla(4);
U_pla1 = output_nompla(5);
eff_pla1 = output_nompla(6);

% Reformating constraints
M(:,4:end) = [];
lb = lb';
ub = ub';
lb(4:end)  = [];
ub(4:end)  = [];
slew_rate = [0.001; 0.1; 0.5];

deltaH = [1e-3 1e-2 1e-2];
Kgrad  = 0.9;
Nsim = 300;
tol  = 1e-5;
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

for i = 1:3
    in_grad = u_ss;
    in_grad(i) = in_grad(i) + deltaH(i);
    
    [x_grad] = OPTIM_SteadyState(in_grad,T_0,T_in,SOFC_data_nominal);
    [~,U_nom2, P_nom2,~,~,~,eff_nom2] = fPrimeMyProject(0,x_grad,in_grad,T_in,SOFC_data_nominal);
    y_ss = [P_nom2; U_nom2];
    [controller, target_config,U_nom3,P_nom3,eff_nom3] = MPC_controller(optimal_var,Ts, x_grad',in_grad,y_ss,Pel_target,T_in,SOFC_data_nominal,M,m,lb,ub,slew_rate);

    x_setpoint = target_config(1:8);
    u_target = target_config(9:11);
    x_plant = initial_config(4:end);
    xhat_plant = x_plant-x_grad';

    u_plant = [];
    P_plant = [];
    U_plant = [];
    eff_plant = [];
    u_previous = initial_config(1:3);

    for j = 1:Nsim
        [u_hat, infeasible] = controller{{xhat_plant(:,j), x_setpoint, u_target, u_previous}};
        u_plant(:,j) = u_hat(:,1)+in_grad;
        u_plant(:,j) = MPC_apply_slewrate(u_plant(:,j),u_previous,slew_rate);
        u_previous = u_plant(:,j);

        sol = ode15s(@(t,x) fPrimeMyProject(t,x,u_plant(:,j),T_in,SOFC_data_plant), [0 Ts], x_plant(:,j), opts);
        x_plant(:,j+1) = sol.y(:,end);
        xhat_plant(:,j+1) = x_plant(:,j+1) - x_grad';
        [dx_states,U_plant(j),P_plant(j),~,~,~,eff_plant(j)] = fPrimeMyProject(0,x_plant(:,j+1),u_plant(:,j),T_in,SOFC_data_plant);

        check_steady_plant = norm(dx_states);
        if check_steady_plant < tol
            break;
        end
    end
    
    P_pla2 = P_plant(end);
    U_pla2 = U_plant(end);
    eff_pla2 = eff_plant(end);
    
    grad_Ucell_nomnl(i) = (U_nom2-U_nom1)/deltaH(i);
    grad_Ucell_plant(i) = (U_pla2-U_pla1)/deltaH(i);
    grad_Power_nomnl(i) = (P_nom2-P_nom1)/deltaH(i);
    grad_Power_plant(i) = (P_pla2-P_pla1)/deltaH(i);
    grad_Effic_nomnl(i) = (eff_nom2-eff_nom1)/deltaH(i);
    grad_Effic_plant(i) = (eff_pla2-eff_pla1)/deltaH(i);
    
    P_plant_hist = [P_plant_hist, P_plant];
    eff_plant_hist = [eff_plant_hist, eff_plant];
    U_plant_hist = [U_plant_hist, U_plant];
    input_plant  = [input_plant u_plant];
    output_plant = [output_plant, x_plant(:,1:end-1)];
    last_config  = [u_plant(:,end); x_plant(:,end-1)];
end

modif(2,1:3) = (1-Kgrad)*modif(2,1:3) + Kgrad*(grad_Ucell_plant-grad_Ucell_nomnl); % Gradient cell voltage
modif(3,1:3) = (1-Kgrad)*modif(3,1:3) + Kgrad*(grad_Power_plant-grad_Power_nomnl); % Gradient power
modif(1,1:3) = (1-Kgrad)*modif(1,1:3) + Kgrad*(grad_Effic_plant-grad_Effic_nomnl); % Gradient efficiency

end