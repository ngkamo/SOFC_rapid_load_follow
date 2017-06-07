function [dxH, outflow] = HEX_Dynamics(t,xH,inflow,T_HEX,SOFC_data)
global myProblem
% inflowh = [H2       H2O       CH4       CO        CO2       O2        N2];
inflowh = [inflow(3) inflow(1) inflow(2) inflow(5) inflow(4) inflow(6) inflow(7)];    
inflowc = [inflow(8) inflow(8)*SOFC_data.CH.r_N2];

Th_in   = T_HEX(1);    % burner temperature 
Tc_in   = T_HEX(2);    % air stream temperature

outflow = inflowh; % outflow of the system
%% HEX parameters
T_ref = SOFC_data.cst.T_ref;
patm  = SOFC_data.cst.p_atm;

% Dimensions
L            = SOFC_data.hex.L;
W            = SOFC_data.hex.W;
H            = SOFC_data.hex.H;
e            = SOFC_data.hex.e;
plate_number = SOFC_data.hex.nPlate;

A        = L*W;
Atot     = A*(plate_number+2);
Vchannel = A*H*plate_number;
Vmetal   = (plate_number-2)*A*e;

% Metal carateristics
CpMet  = SOFC_data.hex.CpMet;
rhoMet = SOFC_data.hex.rhoMet;

%% Data extraction
Th_out = xH(1);
Tc_out = xH(2);

% Fume side
nCH4 = inflowh(1);
nH2  = inflowh(2);
nH2O = inflowh(3);
nCO2 = inflowh(4);
nCO  = inflowh(5);
nO2h  = inflowh(6);
nN2h  = inflowh(7);

% Air side
nO2  = inflowc(1);
nN2  = inflowc(2);

%% Flow properties computations

% Molar fractions
nhtot = nCH4+nH2+nH2O+nCO2+nCO+nO2h+nN2h;
nctot = nO2+nN2;

xh = [nCH4 nH2 nH2O nCO2 nCO nO2h nN2h]/nhtot;
xc = [nO2 nN2]/nctot;

% Average molar masses
Mmh = [SOFC_data.cst.MCH4, SOFC_data.cst.MH2, SOFC_data.cst.MH2O, SOFC_data.cst.MCO2, SOFC_data.cst.MCO, SOFC_data.cst.MO2, SOFC_data.cst.MN2];
Mmc = [SOFC_data.cst.MO2, SOFC_data.cst.MN2];

Mmh = sum(Mmh.*xh);
Mmc = sum(Mmc.*xc);

% Average densities
rhoh = [COR_Density(Th_out,patm,'CH4'),COR_Density(Th_out,patm,'H2'),...
        COR_Density(Th_out,patm,'H2O'),COR_Density(Th_out,patm,'CO2'),...
        COR_Density(Th_out,patm,'CO'),COR_Density(Th_out,patm,'O2'),...
        COR_Density(Th_out,patm,'N2')];
rhoc = [COR_Density(Tc_out,patm,'O2'),COR_Density(Tc_out,patm,'N2')];

rhoMh = sum(rhoh.*xh);
rhoMc = sum(rhoc.*xc);

% Average viscosity
visch = [COR_Viscosity(Th_out,'CH4'),COR_Viscosity(Th_out,'H2'),...
         COR_Viscosity(Th_out,'H2O'),COR_Viscosity(Th_out,'CO2'),...
         COR_Viscosity(Th_out,'CO'),COR_Viscosity(Th_out,'O2'),...
         COR_Viscosity(Th_out,'N2')];
viscc = [COR_Viscosity(Tc_out,'O2'),COR_Viscosity(Tc_out,'N2')];

viscMh = sum(visch.*xh);
viscMc = sum(viscc.*xc);

% Average specific heat 
Cph = [COR_Cp(Th_out,'CH4'),COR_Cp(Th_out,'H2'),...
       COR_Cp(Th_out,'H2O'),COR_Cp(Th_out,'CO2'),...
       COR_Cp(Th_out,'CO'),COR_Cp(Th_out,'O2'),...
       COR_Cp(Th_out,'N2')];
Cpc = [COR_Cp(Tc_out,'O2'),COR_Cp(Tc_out,'N2')];

CpMh = sum(Cph.*xh);
CpMc = sum(Cpc.*xc);

% Average conductivity
kh = [COR_Conductivity(Th_out,'CH4'),COR_Conductivity(Th_out,'H2'),...
      COR_Conductivity(Th_out,'H2O'),COR_Conductivity(Th_out,'CO2'),...
      COR_Conductivity(Th_out,'CO'),COR_Conductivity(Th_out,'O2'),...
      COR_Conductivity(Th_out,'N2')];
kc = [COR_Conductivity(Tc_out,'O2'),COR_Conductivity(Tc_out,'N2')];

kMh = sum(kh.*xh);
kMc = sum(kc.*xc);

%% Heat exchange computation
% Convective heat exchange
Qh = nhtot/(Mmh*rhoMh);
Qc = nctot/(Mmc*rhoMc);

vh = Qh/A; 
vc = Qc/A;

ReMh = L*vh*rhoMh/viscMh;
ReMc = L*vc*rhoMc/viscMc;

Prh = sum((Cph./Mmh).*xh)*viscMh/kMh;
Prc = sum((Cpc./Mmc).*xc)*viscMc/kMc;

if ReMh > 400000
    NuMh = 0.037*nthroot(ReMh,5/4)*nthroot(Prh,3);
elseif ReMh <= 400000 && ReMh >= 0
    NuMh = 0.664*ReMh^(1/2)*nthroot(Prh,3);
else 
    NuMh = 1e10*ReMh*nthroot(Prh,3); 
end

if ReMc > 400000
    NuMc = 0.037*nthroot(ReMc,5/4)*nthroot(Prc,3);
elseif ReMc <= 400000 && ReMc >= 0   
    NuMc = 0.664*ReMc^(1/2)*nthroot(Prc,3);
else
    NuMc = 1e10*ReMc*nthroot(Prc,3);
end

hMh = NuMh*kMh/L;
hMc = NuMc*kMc/L;

Rconvh = 1/(Atot*hMh);
Rconvc = 1/(Atot*hMc);

% Conductive heat exchange (stainless steel)
Ts = (xH(1)+xH(2))/2;
if Ts >= 0
    ks = 1.22*Ts^0.437;
else
    ks = 1e10*Ts;
end
Rcond = e/(ks*Atot*2);

% Total thermal resistance
htot = 1/((Rconvh+Rconvc+Rcond)*Atot);

% Heat tranfert between channels
Q_t = -1*htot*Atot*(xH(1)-xH(2));

%% Inlet and oultet enthalpy flows
Hh  = [SOFC_data.cst.H_CH4, SOFC_data.cst.H_H2, SOFC_data.cst.H_H2O, SOFC_data.cst.H_CO2, SOFC_data.cst.H_CO,SOFC_data.cst.H_O2, SOFC_data.cst.H_N2];
Hc  = [SOFC_data.cst.H_O2, SOFC_data.cst.H_N2];

Dhinh = [COR_Energies(T_ref,Th_in,'CH4'),COR_Energies(T_ref,Th_in,'H2'),...
          COR_Energies(T_ref,Th_in,'H2O'),COR_Energies(T_ref,Th_in,'CO2'),...
          COR_Energies(T_ref,Th_in,'CO'),COR_Energies(T_ref,Th_in,'O2'),...
          COR_Energies(T_ref,Th_in,'N2')];    
nHinh = (Dhinh+Hh).*[nCH4 nH2 nH2O nCO2 nCO nO2h nN2h];

Dhinc  = [COR_Energies(T_ref,Tc_in,'O2'),COR_Energies(T_ref,Tc_in,'N2')];
nHinc  = (Dhinc+Hc).*[nO2 nN2];

Dhouth = [COR_Energies(T_ref,Th_out,'CH4'),COR_Energies(T_ref,Th_out,'H2'),...
          COR_Energies(T_ref,Th_out,'H2O'),COR_Energies(T_ref,Th_out,'CO2'),...
          COR_Energies(T_ref,Th_out,'CO'),COR_Energies(T_ref,Th_out,'O2'),...
          COR_Energies(T_ref,Th_out,'N2')];
nHouth = (Dhouth+Hh).*[nCH4 nH2 nH2O nCO2 nCO nO2h nN2h];

Dhoutc = [COR_Energies(T_ref,Tc_out,'O2'),COR_Energies(T_ref,Tc_out,'N2')];
nHoutc = (Dhoutc+Hc).*[nO2 nN2];

hinh  = sum(nHinh);
hinc  = sum(nHinc);
houth = sum(nHouth);
houtc = sum(nHoutc);

%% Energy balance
dxH    = zeros(2,1);

if myProblem.TC.TimeConstantOff == 0
    dxH(1) = (1/100)*1/(Vchannel*CpMh*rhoMh/Mmh)*(Q_t + hinh-houth);
    dxH(2) = (1/100)*1/(Vchannel*CpMc*rhoMc/Mmc)*(-Q_t + hinc-houtc);
else
    dxH(1) = (Q_t + hinh-houth);
    dxH(2) = (-Q_t + hinc-houtc);
end

end
