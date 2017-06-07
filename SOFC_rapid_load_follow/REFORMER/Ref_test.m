clear all
close all
clc

global SOFC_data_nominal myProblem


N_c = 60;
prjname_SOFC       = 'Solid Oxide Fuel Cell';
SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);

prjname_myProblem  = 'Real Time Optimization';
myProblem = data_myProblem(prjname_myProblem);

% T_f = 298.15;
% SOFC_data_nominal.CE.T_f = T_f;

% nCH4in = inflow(1);
% nH2Oin = inflow(2);
inflow = 0.0042*[1 2.5];
% inflow = 0.00650833496652208*[1 2.1];
% inflow = 0.00475494135060608*[1 2.1];

% TCH4in = Tflow(1);
% TH2Oin = Tflow(2);

Tflow = [473.15 473.15];

T_burn = 1150+273.15;
 
x_0 = 250+273.15;
options = optimset('Algorithm','trust-region-dogleg',...
                                    'LargeScale','on',...
                                    'Diagnostics','on',...
                                    'Display','on' ,...
                                    'FunValCheck','on',...
                                    'MaxFunEvals',10000,...
                                    'MaxIter',10000,...
                                    'TolFun',1e-9,...
                                    'TolX',1e-9);

myProblem.TC.TimeConstantOff = 1;

Tr = fsolve(@(x) ref_rate_Tr(0,x,inflow,Tflow,T_burn,SOFC_data_nominal),x_0,options);
Tr1 = Tr - 273.15
myProblem.TC.TimeConstantOff = 0;







TSPAN = 0:10:2000;
Y0 = x_0;
tic
[TOUT,YOUT] = ode23s(@(t,x) ref_rate_Tr(t,x,inflow,Tflow,T_burn,SOFC_data_nominal),TSPAN,Y0);
toc

[n_y,~] = size(YOUT);
for i = 1:1:n_y
[~,outflow(i,:),~,~,~] = ref_rate_Tr(TOUT(i),YOUT(i,:),inflow,Tflow,T_burn,SOFC_data_nominal);
end
%% Plotting
close all
% load('REF.mat')

figure
plot(TOUT/60,YOUT-273.15,'LineWidth',2); hold on
xlabel('time [min]','FontSize',15)
ylabel('Reformer temperature [{\circ}C]','FontSize',15)
grid on

nH2out  = outflow(1);
nH2Oout = outflow(2);
nCH4out = outflow(3);
nCOout  = outflow(4);
nCO2out = outflow(5);



% figure
% subplot(2,1,1)
% plot(TOUT/60,outflow(:,1),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' concentration [mol/l] ','FontSize',15)
% legend({'nCH4out'},'FontSize',15)
% grid on
% subplot(2,1,2)
% plot(TOUT/60,outflow(:,3),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' concentration [mol/l] ','FontSize',15)
% legend({'nH2Oout'},'FontSize',15)
% grid on
% 
% figure
% subplot(3,1,1)
% plot(TOUT/60,outflow(:,2),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' concentration [mol/l] ','FontSize',15)
% legend({'nH2out'},'FontSize',15)
% grid on
% subplot(3,1,2)
% plot(TOUT/60,outflow(:,5),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' concentration [mol/l] ','FontSize',15)
% legend({'nCOout'},'FontSize',15)
% grid on
% subplot(3,1,3)
% plot(TOUT/60,outflow(:,4),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' concentration [mol/l] ','FontSize',15)
% legend({'CO2out'},'FontSize',15)
% grid on

% figure
% plot(TOUT/60,outflow(:,3)./outflow(:,1),'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' steam to carbon ratio [-] ','FontSize',15)
% legend({'S/C'},'FontSize',15)
% grid on

% ntot = repmat(sum(outflow')',1,5);
% figure
% plot(TOUT/60,outflow./ntot,'LineWidth',2)
% xlabel(' time [min] ','FontSize',15)
% ylabel(' composition [-] ','FontSize',15)
% legend({'nCH4out','nH2out','nH2Oout','nCO2out','nCOout'},'FontSize',15)
% grid on

figure
plot(TOUT/60,100*(inflow(1)-outflow(:,3))./inflow(1),'LineWidth',2)
xlabel(' time [min] ','FontSize',15)
ylabel(' conversion [%] ','FontSize',15)
grid on

% figure
% plot(YOUT-273.15,100*(inflow(1)-outflow(:,1))./inflow(1),'LineWidth',2);
% hold on
% plot([600 650 750],[55 70 100],'ro','LineWidth',2)
% ylabel(' conversion [%] ','FontSize',15)
% xlabel(' Reformer temperature [{\circ}C] ','FontSize',15)
% grid on
