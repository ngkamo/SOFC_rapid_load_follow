function [dx,isterm,dir] = event(~,x,in_opt,T_in,SOFC_data_plant)
    [dx_global,~,~,~,~,~,~] = fPrimeMyProject(0,x,in_opt,T_in,SOFC_data_plant);
    dx = norm(dx_global) - 1e-4;
    isterm = 1;
    dir = 0;
end