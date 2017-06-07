function [inflow_A,outflow,react,r_Ref] = SOFC_MassBalance(inflow,u,SOFC_data)
N_c   = SOFC_data.CE.N_c;
F     = SOFC_data.cst.F;
I     = u(3);

nH2_outREF  = inflow(1);
nH2O_outREF = inflow(2);
nCH4_outREF = inflow(3);
nCO_outREF  = inflow(4);
nCO2_outREF = inflow(5);
nO2_outHEX  = inflow(6);
nN2_outHEX  = inflow(7);

r_Ref = nCH4_outREF;

nH2in = nH2_outREF + 4*r_Ref;
nH2Oin = nH2O_outREF - 2*r_Ref;
nCH4in = nCH4_outREF - r_Ref;
nCOin  = nCO_outREF;
nCO2in = nCO2_outREF + r_Ref;

nO2in  = nO2_outHEX;
nN2in  = nN2_outHEX;

nH2r  = (N_c*I)/(2*F);

nH2Or = (N_c*I)/(2*F);
nCH4r = 0;
nCOr  = 0;
nCO2r = 0;
nO2r  = (N_c*I)/(4*F);
if nO2r > nO2in
    nO2r = nH2in/2;
end
nN2r  = 0;

% Output 
nH2  = [ nH2in ,  nH2r,  nH2in -  nH2r]; 
nH2O = [nH2Oin , nH2Or, nH2Oin + nH2Or];
nCH4 = [nCH4in , nCH4r, nCH4in - nCH4r];
nCO  = [ nCOin ,  nCOr,  nCOin -  nCOr];
nCO2 = [nCO2in , nCO2r, nCO2in - nCO2r];
nO2  = [ nO2in ,  nO2r,  nO2in -  nO2r];
nN2  = [ nN2in ,  nN2r,  nN2in -  nN2r];


inflow_A = [nH2(1), nH2O(1), nCH4(1), nCO(1), nCO2(1), nO2(1), nN2(1)];
react    = [nH2(2), nH2O(2), nCH4(2), nCO(2), nCO2(2), nO2(2), nN2(2)];
outflow  = [nH2(3), nH2O(3), nCH4(3), nCO(3), nCO2(3), nO2(3), nN2(3)];
end


