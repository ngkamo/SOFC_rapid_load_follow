function [dxB,outflow] = BURN_Dynamics(t, xB,Tf_in, Ta_in, inflow,Tref,SOFC_data)
global myProblem 

%% Burner parameters
R       = SOFC_data.cst.R;
T_ref   = SOFC_data.cst.T_ref;
Vburner = SOFC_data.burn.V;
patm    = SOFC_data.cst.p_atm;

%% Data extraction
Tb = xB;

% Fuel side
nCH4in = inflow(3);
nH2in  = inflow(1);
nH2Oin = inflow(2);
nCO2in = inflow(5);
nCOin  = inflow(4);

% Air side
nO2in  = inflow(6);
nN2in  = inflow(7);

%% Mass balance
nN2out  = nN2in;
nCO2out = nCO2in+nCH4in+nCOin;
nH2Oout = nH2Oin+2*nCH4in+nH2in;
nO2out  = nO2in-2*nCH4in-1/2*nCOin-1/2*nH2in;

% outflow = [0 0 nH2Oout nCO2out 0 nO2out nN2out];
%outflow = [H2  H2O  CH4 CO  CO2     O2     N2];
outflow =  [0 nH2Oout 0  0  nCO2out nO2out nN2out];

%% Flow properties computations
% Molar fractions
ntotout = nH2Oout+nCO2out+nO2out+nN2out;
xout = [nH2Oout nCO2out nO2out nN2out]/ntotout;

% Mass inside the burner
ntot = patm*Vburner/(Tb*R);

% Average specific heat 
Cp = [COR_Cp(Tb,'H2O'),COR_Cp(Tb,'CO2'),...
      COR_Cp(Tb,'O2'),COR_Cp(Tb,'N2')];
CpM = sum(Cp.*xout);

%% Inlet and oultet enthalpy flows
% Inlet fuel side
HinF  = [SOFC_data.cst.H_CH4,SOFC_data.cst.H_H2,SOFC_data.cst.H_H2O,SOFC_data.cst.H_CO2,SOFC_data.cst.H_CO];
DhinF = [COR_Energies(T_ref,Tf_in,'CH4'),COR_Energies(T_ref,Tf_in,'H2'),...
        COR_Energies(T_ref,Tf_in,'H2O'),COR_Energies(T_ref,Tf_in,'CO2'),...
        COR_Energies(T_ref,Tf_in,'CO')];    
nHinF = (DhinF+HinF).*[nCH4in nH2in nH2Oin nCO2in nCOin];

% Inlet air side
HinA   = [SOFC_data.cst.H_O2 SOFC_data.cst.H_N2];
DhinA  = [COR_Energies(T_ref,Ta_in,'O2'),COR_Energies(T_ref,Ta_in,'N2')];
nHinA  = (DhinA+HinA).*[nO2in nN2in];

% Outlet burner
HoutB  = [SOFC_data.cst.H_O2,SOFC_data.cst.H_N2,SOFC_data.cst.H_CO2,SOFC_data.cst.H_H2O];
DhoutB = [COR_Energies(T_ref,Tb,'O2'),COR_Energies(T_ref,Tb,'N2'),...
         COR_Energies(T_ref,Tb,'CO2'),COR_Energies(T_ref,Tb,'H2O')];
nHoutB = (DhoutB+HoutB).*[nO2out nN2out nCO2out nH2Oout];

hinF  = sum(nHinF);
hinA  = sum(nHinA);
houtB = sum(nHoutB);

DMet = -0.8026*1e9/1e3;
Dhy  = -0.2418*1e9/1e3;
DhyCO = -0.2830*1e9/1e3;

Dtot = nCH4in*DMet + nH2in*Dhy + nCOin*DhyCO;

A_r = 0.001; % m^2
hg = 1000/10; % W/(K*m^2)
Q_convec = hg*A_r*(Tb - Tref);

%% Energy balance

if myProblem.TC.TimeConstantOff == 0    
    dxB = (1/1)*1/(ntot*CpM)*(hinA+hinF-houtB - Dtot - Q_convec);    
else
    dxB = hinA+hinF-houtB - Dtot - Q_convec;
end

end