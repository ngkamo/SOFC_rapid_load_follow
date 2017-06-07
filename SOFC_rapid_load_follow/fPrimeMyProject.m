function [dx_global,Ucell,Pel,FU,L_air,eta_el,eta_sys] = fPrimeMyProject(t,x_global,u_s,T_in,SOFC_data)
x_REF  = x_global(1);
x_SOFC = x_global(2:5);
x_BURN = x_global(6);
x_HEX  = x_global(7:8);

%% converting L/min into mol/s
R = 8.314462175;    % [J/K/mol] 
P = 100e3;          % [Pa]
T_st = 273.15;      % [K]

V_O2cathinp = u_s(2)*0.21;

nCH4inp    = (P*u_s(1))/((R*T_st)*(6e+4));      % [mol/s]     methane flow rate
nO2cathinp = (P*V_O2cathinp)/((R*T_st)*(6e+4)); % [mol/s]     oxygen flow rate

input = [nCH4inp nO2cathinp u_s(3)];

%%
TCH4in = T_in(1);
TH2Oin = T_in(2); 

nCH4in    = input(1);
nH2Oin    = 2.5*input(1);

%% Reformer
[dx_REF,outflow,~,~,~] = ref_rate_Tr(t,x_REF,[nCH4in,nH2Oin],[TCH4in,TH2Oin],x_BURN,SOFC_data);

%% SOFC
[dx_SOFC,Ucell,Pel,FU,L_air,eta_el,outflow_SOFC] = SOFC_Dynamics(t,x_SOFC,input,outflow,[x_REF x_HEX(2)],SOFC_data);
                                                   
eta_sys = Pel/(SOFC_data.cst.LHVCH4) * (1 / (SOFC_data.cst.MCH4 *input(1)));
%% Burner
[dx_BURN,outflowb] = BURN_Dynamics(t,x_BURN,x_SOFC(3),x_SOFC(4),outflow_SOFC ,x_REF,SOFC_data);

%% Heat exchanger air inlet
[dx_HEX, ~] = HEX_Dynamics(t,x_HEX,[outflowb,input(2)],[x_BURN T_in(3)],SOFC_data);

%% Global states
dx_global        = zeros(8,1);
dx_global(1)     = (1/1)*dx_REF;
dx_global(2:5)   = (1/1)*dx_SOFC;
dx_global(6)     = (1/1)*dx_BURN;
dx_global(7:8)   = (1/1)*dx_HEX';
end