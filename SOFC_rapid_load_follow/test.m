close all
timestep = Ts*[0:Nsim]/60;
set(0,'DefaultFigureWindowStyle','docked')
figure('Name','power')
[ts,ys] = stairs(timestep(1:end-1),y(1,:)');
plot(ts,ys)
xlabel('Time [min]')
ylabel('Power [W]')
title('Power [W]')
legend('linear','nonlinear')
grid on


figure('Name','voltage')
[ts,ys] = stairs(timestep(1:end-1),y(2,:));
plot(ts,ys)
xlabel('Time [min]')
ylabel('Voltage [V]')
title('Voltage [V]')
legend('linear')
grid on

% figure('Name','efficiency')
% [ts,ys] = stairs(timestep(1:end-1),y(3,:));
% plot(ts,ys, timestep(1:end-1),eta_nonlin);
% xlabel('Time [min]')
% ylabel('Efficiency [-]')
% title('Efficiency')
% legend('linear','nonlinear')
% grid on

figure('Name','methane')
[ts,ys] = stairs(timestep(1:end-1),u_lin(1,:)+u_opt(target,1));
plot(ts,ys)
xlabel('Time [min]')
ylabel('q_{CH4} [L/min]')
title('Input methane [L/min]')
legend('linear')
grid on

figure('Name','air')
[ts,ys] = stairs(timestep(1:end-1),u_lin(2,:)+u_opt(target,2));
plot(ts,ys)
xlabel('Time [min]')
ylabel('q_{air} [L/min]')
title('Input air flow [L/min]')
legend('linear')
grid on

figure('Name','current')
[ts,ys] = stairs(timestep(1:end-1),u_lin(3,:)+u_opt(target,3));
plot(ts,ys);
xlabel('Time [min]')
ylabel('current I [A]')
title('Current [A]')
legend('linear')
grid on

figure('Name','temperatures')
[ts,ys] = stairs(timestep(1:end),x_hat(3,:));
plot(ts,ys);
xlabel('Time [min]')
ylabel('current I [A]')
title('Current [A]')
legend('linear')
grid on