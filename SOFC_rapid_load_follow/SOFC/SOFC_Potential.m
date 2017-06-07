function [U_Nernst,n_loss] = SOFC_Potential(xt,u,outflow,T_a,FU,L_air,SOFC_data)
R          = SOFC_data.cst.R;
F          = SOFC_data.cst.F;

k0         = SOFC_data.kinetic.k0;
E_actcath  = SOFC_data.kinetic.E_actcath;
E_disscath = SOFC_data.kinetic.E_disscath;
sig0_el    = SOFC_data.kinetic.sig0_el;
E_el       = SOFC_data.kinetic.E_el;
R0         = SOFC_data.kinetic.R0;
A_act      = SOFC_data.CE.A_act;
h_el       = SOFC_data.EL.h;
R_MIC      = SOFC_data.IC.R_MIC;
T_fuel     = T_a(1); % at the inlet
FU_adj     = SOFC_data.CE.FU_adj;

I    = u(3);

x_elchem = .7479; 
T_elchem = T_fuel + (xt(3) - T_fuel)*x_elchem;

nH2  = outflow(1);
nH2O = outflow(2);
nO2  = outflow(6);
nN2  = outflow(7);

xH2anout   = nH2/(nH2 + nH2O);
xH2Oanout  = nH2O/(nH2 + nH2O);
xO2cathout = nO2/(nO2 + nN2) ;
xN2cathout = nN2/(nO2 + nN2);

x_out = [xH2anout xH2Oanout xO2cathout xN2cathout];

h_out = COR_Enthalpy(xt(3),xt(4));
s_out = COR_Entropy(xt(3),xt(4),x_out);
g_out = h_out - xt(3)*s_out;

G_out = g_out(2) - g_out(1) - 0.5*g_out(3); 
G_in  = 0; 
dG    = G_out - G_in;

U_Nernst = -dG/2/F;

% Overpotential losses

% Activation overpotential :
i         = I/A_act; 
i0_cath   = (R*T_elchem*2/F)*k0*exp(-E_actcath/R/T_elchem); 
n_actcath = R*T_elchem/F*asinh(i/(2*i0_cath));

% Dissociation overpotential :
n_diss = i*R0*(0.21)^(-0.5)*exp(E_disscath/(R*T_elchem)); 

% Ionic overpotential
sig_el = sig0_el*exp(-E_el/R/T_elchem);
n_el = i*(h_el/sig_el);

% Concentration overpotential
n_diffan   = -R*T_elchem/2/F*log(1-(FU+FU_adj));

n_diffcath = -R*T_elchem/2/F*log(1- FU/L_air);

% Ohmic overpotential
MIC_loss = R_MIC*i*1e-4; 

% Overpoential losses vector output
n_loss = [n_actcath, n_diss, n_el, n_diffan, n_diffcath, MIC_loss];

end

