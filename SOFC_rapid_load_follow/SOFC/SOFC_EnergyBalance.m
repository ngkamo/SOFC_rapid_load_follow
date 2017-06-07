function [dT_electrolyte, dT_interconnector, dT_fuel, dT_air] = SOFC_EnergyBalance(xt,u,Ucell,inflow_A,outflow,react,T_in,r_Ref,SOFC_data)
global myProblem
sigma  = SOFC_data.cst.sigma;
e_VF   = SOFC_data.CE.e_VF;
A_loss = SOFC_data.CE.A_loss;
R      = SOFC_data.cst.R;
p_atm  = SOFC_data.cst.p_atm;
N_c    = SOFC_data.CE.N_c;

nH2_out  = outflow(1);
nO2_out  = outflow(6);
nN2_out  = outflow(7);
nH2O_out = outflow(2);
nCH4_out = 0;%
nCO_out  = 0;%
nCO2_out = 0;%

nH2_reac  = react(1);
nO2_reac  = react(6);
nH2O_reac = react(2);

nH2_in  = inflow_A(1);
nO2_in  = inflow_A(6);
nN2_in  = inflow_A(7);
nH2O_in = inflow_A(2);
nCH4_in = 0;%
nCO_in  = 0;% 
nCO2_in = 0;%

T_ref   = SOFC_data.cst.T_ref;

Tcin    = T_in(2);
Tfin    = T_in(1);

% Electrode data
h_el   = (SOFC_data.AN.h+SOFC_data.CA.h+SOFC_data.EL.h)*N_c ;  % [m]
A_el   = SOFC_data.AN.l_x*SOFC_data.AN.l_y;  % [m^2] Anode, cathode and electrolyte are assumed to have the same area
rho_el = (SOFC_data.AN.h*SOFC_data.AN.rho+SOFC_data.CA.h*SOFC_data.CA.rho+SOFC_data.EL.h*SOFC_data.EL.rho)/(SOFC_data.AN.h+SOFC_data.CA.h+SOFC_data.EL.h); % [kg/m3] Average density for anode, cathode and electrolyte
cp_el  = (SOFC_data.AN.C_p*SOFC_data.AN.h+SOFC_data.CA.h*SOFC_data.CA.C_p+SOFC_data.EL.h*SOFC_data.EL.C_p)/(SOFC_data.AN.h+SOFC_data.CA.h+SOFC_data.EL.h); % [J/kg/K] Average Cp for anode, cathode and electrolyte
e_a    = SOFC_data.AN.e_a;   % [-] 
e_c    = SOFC_data.CA.e_c;   % [-] 

% Interconnector data
h_i   = SOFC_data.IC.h*(N_c+1); % [m]
A_i   = SOFC_data.IC.l_x*SOFC_data.IC.l_y;  % [m^2]
rho_i = SOFC_data.IC.rho;  % [kg/m3]
cp_i  = SOFC_data.IC.C_p;  % [J/kg/K]
e_i   = SOFC_data.IC.e_i;  % [-]

% Fuel channel data
V_channel = SOFC_data.CH.V_a*N_c;

% Convection coefficients
na_in = nH2_in+nH2O_in+nCH4_in+nCO_in+nCO2_in;
Vadot = na_in*R*Tfin/p_atm;
va    = Vadot/(SOFC_data.CH.h_a*SOFC_data.CH.w_a);

nc_in = nO2_in + nN2_in;
Vcdot = nc_in*R*Tcin/p_atm;
vc    = Vcdot/(SOFC_data.CH.h_c*SOFC_data.CH.w_c);

hFe = 5000; %10.45 - va + 10*va^(1/2);
hAe = 5000;%10.45 - vc + 10*vc^(1/2);
hFi = 5000;%10.45 - va + 10*va^(1/2);
hAi = 5000;%10.45 - vc + 10*vc^(1/2);

% Energy Balance Stack
% Electrode balance
% Diffusive heat transfers
QdFin = COR_Energies(T_ref,xt(3),'H2')*nH2_reac;
QdFout = COR_Energies(T_ref,xt(3),'H2O')*nH2O_reac;
QdAin = COR_Energies(T_ref,xt(4),'O2')*nO2_reac;
QdAout = 0;

% Convective heat tranfers
QhFe = (xt(1)-xt(3))*hFe*A_el;
QhAe = (xt(1)-xt(4))*hAe*A_el;

% Radiative losses with air and fuel
QrA = sigma*A_el*(xt(1)^4-xt(2)^4)/(1/e_a+1/e_i-1);
QrF = sigma*A_el*(xt(1)^4-xt(2)^4)/(1/e_c+1/e_i-1);

% Heat of reaction
H_H2O = -241827 + COR_Energies(T_ref,xt(3),'H2O');
H_H2  = 0 + COR_Energies(T_ref,xt(3),'H2');
H_O2  = 0 + COR_Energies(T_ref,xt(3),'O2');
DG    = H_H2O-H_H2-H_O2*1/2;
Qr = nH2_reac*DG;

% Electrical power
W = u(3)*Ucell*N_c;

% Balance equation

if myProblem.TC.TimeConstantOff == 0
    dT_electrolyte = 1/(rho_el*A_el*h_el*cp_el)*((QdFin-QdFout)+(QdAin-QdAout)-(QhFe+QrF+QhAe+QrA)-Qr-W);
else
    dT_electrolyte = (QdFin-QdFout)+(QdAin-QdAout)-(QhFe+QrF+QhAe+QrA)-Qr-W;
end

% Interconnector balance
% Convective heat transfert
QhFi = (xt(3)-xt(2))*hFi*A_i;
QhAi = (xt(4)-xt(2))*hAi*A_i;

% Radiative heat transfert
% Q_loss = e_VF*sigma*(xt(2)^4 - T_f^4)*A_loss;
Q_loss = 0;

% Balance equation
if myProblem.TC.TimeConstantOff == 0
    dT_interconnector = 1/(rho_i*A_i*h_i*cp_i)*(QhFi+QhAi+QrA+QrF-Q_loss);
else
    dT_interconnector = QhFi+QhAi+QrA+QrF-Q_loss;
end

% Fuel side
% reforming reactions
D_RefShif = 164.95e3;
Q_refor = r_Ref*D_RefShif;

% Heat in
QFin = nH2_in *COR_Energies(T_ref,Tfin,'H2') +nCO_in *COR_Energies(T_ref,Tfin,'CO') +...
       nCO2_in*COR_Energies(T_ref,Tfin,'CO2')+nCH4_in*COR_Energies(T_ref,Tfin,'CH4')+...
       nH2O_in*COR_Energies(T_ref,Tfin,'H2O');
   
% Heat out
QFout = nH2_out *COR_Energies(T_ref,xt(3),'H2') +nCO_out *COR_Energies(T_ref,xt(3),'CO') +...
        nCO2_out*COR_Energies(T_ref,xt(3),'CO2')+nCH4_out*COR_Energies(T_ref,xt(3),'CH4')+...
        nH2O_out*COR_Energies(T_ref,xt(3),'H2O') ;
    
% Balance equation
na_out = nH2_out+nH2O_out+nCH4_out+nCO_out+nCO2_out;
xa_out = [nH2_out nH2O_out nCH4_out nCO_out nCO2_out]./na_out;
pa_out = xa_out.*p_atm; % Pressure assumed to be constant at atmospheric pressure
cpa_out = [COR_Cp(xt(3),'H2') COR_Cp(xt(3),'H2O') COR_Cp(xt(3),'CH4') COR_Cp(xt(3),'CO') COR_Cp(xt(3),'CO2')];
d = sum(pa_out.*cpa_out);

if myProblem.TC.TimeConstantOff == 0
    dT_fuel = 1/(V_channel/(R*xt(3))*d)*((QFin-QFout)+(QdFout-QdFin)+(QhFe-QhFi)-Q_refor);
else
    dT_fuel = (QFin-QFout)+(QdFout-QdFin)+(QhFe-QhFi)-Q_refor;
end

% Air side
% Heat in
QAin = nO2_in*COR_Energies(T_ref,Tcin,'O2')+nN2_in*COR_Energies(T_ref,Tcin,'N2');

% Heat out
QAout = nO2_out*COR_Energies(T_ref,xt(4),'O2')+nN2_out*COR_Energies(T_ref,xt(4),'N2');

% Balance equation
nc_out = nO2_out + nN2_out;
xc_out = [nO2_out nN2_out]./nc_out;
pc_out = xc_out*p_atm;
cpc_out = [COR_Cp(xt(4),'O2') COR_Cp(xt(4),'N2')];
e = sum(pc_out.*cpc_out);

if myProblem.TC.TimeConstantOff == 0
    dT_air = 1/(V_channel*e/(R*xt(4)))*((QAin-QAout)+(QdAout-QdAin)+(QhAe-QhAi));
else
    dT_air = (QAin-QAout)+(QdAout-QdAin)+(QhAe-QhAi);
end

end

