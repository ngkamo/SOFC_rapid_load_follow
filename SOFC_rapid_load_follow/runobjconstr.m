function [x,f,eflag,outpt,lambda] = runobjconstr(u0,modif,T_in,A,B,lb,ub,SOFC_data,opts)

if nargin == 8
   opts = [];
end

xLast = [];
myf   = [];
myc   = [];
myceq = [];

fun  = @(x) objfun(x,modif,T_in,SOFC_data);
cfun = @(x) constr(x,modif,T_in,SOFC_data);


[x,f,eflag,outpt,lambda] = fmincon(fun,u0,A,B,[],[],lb,ub,cfun,opts);

    function y = objfun(x,modif,T_in,SOFC_data)
        if ~isequal(x,xLast)
            [myf,myc,myceq] = computeall(x,modif,T_in,SOFC_data);
            xLast = x;
        end
        y = myf;
    end


    function [c,ceq] = constr(x,modif,T_in,SOFC_data)
        if ~isequal(x,xLast)
            [myf,myc,myceq] = computeall(x,modif,T_in,SOFC_data);
            xLast = x;
        end
        c = myc;
        ceq = myceq;
    end

end