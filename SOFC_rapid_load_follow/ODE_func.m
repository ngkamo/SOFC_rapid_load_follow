function [dx] = ODE_func(t,x,u,T_in,SOFC_data,K,x_ss,u_ss,var_target)
%     x_hat = x-x_ss;
%     u_hat = -K*(x_hat-var_target(1:8))+var_target(9:11)
%     if u_hat(1)<0
%         u_hat(1) = -u_ss(1);
%     end
%     if u_hat(3)<0
%         u_hat(3) = -u_ss(3);
%     end
%     u = u_hat+u_ss
    [dx,~,~,~,~,~,~] = fPrimeMyProject(t,x,u,T_in,SOFC_data);
end