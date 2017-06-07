function [dx,Ucell,Pel,FU,L_air,eta_el,outflow] = SOFC_Dynamics(t,x,u,inflow,T_in,SOFC_data)

xt = x;

% Temperatures
T_a(1)  = T_in(1);
T_c(1)  = T_in(2);

%% Mass Balance
[inflow_A,outflow,react,r_Ref] = SOFC_MassBalance([inflow, u(2), u(2)*SOFC_data.CH.r_N2],u,SOFC_data);

FU   =   u(3)*SOFC_data.CE.N_c/8/SOFC_data.cst.F/u(1); % Fuel utilisation [-]
L_air =   u(2)/2/u(1); % Lambda air [-]

%% Potential
[U_Nernst,n_loss] = SOFC_Potential(xt,u,outflow,T_a,FU,L_air,SOFC_data);

Ucell = U_Nernst - sum(n_loss);
Pel   = Ucell*SOFC_data.CE.N_c*u(3);

%% Energy Balance
[dT_electrolyte, dT_interconnector, dT_fuel, dT_air] = SOFC_EnergyBalance(xt,u,Ucell,inflow_A,outflow,react,T_in,r_Ref,SOFC_data);

dxt(1) = dT_electrolyte;
dxt(2) = dT_interconnector;
dxt(3) = dT_fuel;
dxt(4) = dT_air;

%% Copy of state
dx = zeros(4,1);
dx = (1/0.10)* dxt';

%% Efficiency
eta_el = (Ucell*SOFC_data.CE.N_c*u(3))/(SOFC_data.cst.LHV) * (1 / (SOFC_data.cst.MH2 * inflow(1)));

end

