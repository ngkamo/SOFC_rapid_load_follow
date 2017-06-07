function [dTr,outflow,rR,ksi,eff] = ref_rate_Tr(t,x,inflow,Tflow,T_burn,SOFC_data)
global myProblem
Tr    = x;
    
TCH4in = Tflow(1);
TH2Oin = Tflow(2);
nCH4in = inflow(1);
nH2Oin = inflow(2);

R     = 8.314; %J/(mol.K)
T_ref = 298.15;
DHr   = 206.1e3;
DHs   = -41.15e3;
P  = 1;

x0 = [0.49 0.99];
[X,FVAL,EXITFLAG] = fsolve(@(x) equi_ref(x,nCH4in,Tr,P),x0,myProblem.Solver.options);                                

rR  = X(1)*nCH4in;
ksi = nCH4in*X(1)*(1-X(2)); 

nCH4out = nCH4in - rR;
nH2Oout = 2.1*nCH4in - rR - ksi;
nCOout  = rR - ksi;
nH2out  = 3*rR + ksi;
nCO2out = ksi;

outflow = [nH2out nH2Oout nCH4out nCOout nCO2out];

A_r = 0.001; % m^2
Q_rad = 0;
hg = 1000/10; % W/(K*m^2)
Q_convec = hg*A_r*(Tr - T_burn);

dTr   = nCH4in*energy(TCH4in,T_ref,1) + nH2Oin*energy(TH2Oin,T_ref,2) - ...
       (nCH4out*energy(Tr,T_ref,1) + nH2Oout*energy(Tr,T_ref,2) + nH2out*energy(Tr,T_ref,3) + nCO2out*energy(Tr,T_ref,4) + nCOout*energy(Tr,T_ref,5)) ... 
        - rR*DHr - ksi*DHs - Q_rad - Q_convec;

if myProblem.TC.TimeConstantOff == 0
    dTr = (1/0.1)*1/(SOFC_data.ref.rho*SOFC_data.ref.Cp*SOFC_data.ref.Vr)*dTr; 
else
    dTr = dTr;
end

HHV_H2  = 286; % kJ/mol
HHV_CH4 = 889; % kJ/mol
eff = nH2out*HHV_H2/((nCH4in+rR)*HHV_CH4);
end
