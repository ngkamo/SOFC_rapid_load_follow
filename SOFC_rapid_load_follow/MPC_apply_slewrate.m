function inp_nonlin = MPC_apply_slewrate(inp_nonlin,inp_previous,slew_rate)
 if inp_nonlin(1)-inp_previous(1) > slew_rate(1)
    inp_nonlin(1) = inp_previous(1) + slew_rate(1);
elseif inp_nonlin(1)-inp_previous(1) < -slew_rate(1)
    inp_nonlin(1) = inp_previous(1) - slew_rate(1);
end
if (inp_nonlin(2)-inp_previous(2)) > slew_rate(2)
    inp_nonlin(2) = inp_previous(2) + slew_rate(2);
elseif inp_nonlin(2)-inp_previous(2) < -slew_rate(2)
    inp_nonlin(2) = inp_previous(2) - slew_rate(2);
end
if inp_nonlin(3)-inp_previous(3) > slew_rate(3)
    inp_nonlin(3) = inp_previous(3) + slew_rate(3);
elseif inp_nonlin(3)-inp_previous(3) < -slew_rate(3)
    inp_nonlin(3) = inp_previous(3) - slew_rate(3);
end
end