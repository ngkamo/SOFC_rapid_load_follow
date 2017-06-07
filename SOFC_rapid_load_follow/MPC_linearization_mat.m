function [A,B,C,D] = MPC_linearization_mat(u_ss,dx_ss,Pel_ss,U_ss,T_in,SOFC_data)
deltaA = 0.001*ones(1,8);
differentialA = [];
component = 8;
A = zeros(8,8);
B = zeros(8,3);

for i = 1:length(deltaA)
    u_grad = u_ss;
    u_grad(3+i) = u_grad(3+i) + deltaA(i);
    
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data);
    differentialA(:,1) = (dx_test-dx_ss)/deltaA(i);
    
    A(:,i) = differentialA;
end

deltaB = [1e-6 1e-5 1e-5];
differentialB = [];
for i = 1:length(deltaB)
    u_grad = u_ss;
    u_grad(i) = u_grad(i) + deltaB(i);
    
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data);
    differentialB(:,1) = (dx_test-dx_ss)/deltaB(i);
    
    B(:,i) = differentialB;
end

% C & D
deltaCPel = 1e-2*ones(1,8);
deltaCUc  = 1e-2*ones(1,8);
deltaCEff = 1e-4*ones(1,8);
differentialCPel = [];
differentialCUc  = [];
differentialCEff = [];
C = zeros(2,8);

for i = 1:length(deltaCPel)
    u_grad = u_ss;
    u_grad(3+i) = u_grad(3+i) + deltaCPel(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data);
    C(1,i) = (Pel_test-Pel_ss)/deltaCPel(i);
    C(2,i) = (Ucell_test-U_ss)/deltaCUc(i);
end

deltaDPel = [1e-8 1e-6 1e-6];
deltaDUc  = [1e-8 1e-2 1e-3];
deltaDEff = [1e-7 1e-2 1e-2];
differentialD = [];
D = zeros(2,3);
for i = 1:length(deltaDPel)
    u_grad = u_ss;
    u_grad(i) = u_grad(i) + deltaDPel(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data);
    D(1,i) = (Pel_test-Pel_ss)/deltaDPel(i);
    
    u_grad = u_ss;
    u_grad(i) = u_grad(i) + deltaDUc(i);
    [dx_test,Ucell_test,Pel_test,FU_test,Lair_test,~,~] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data);
    D(2,i) = (Ucell_test-U_ss)/deltaDUc(i);
end

end