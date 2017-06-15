clear all;
close all;
clc;
warning on
addpath('REFORMER','SOFC','COR','HEX','BURN');

%%
global myProblem Ps_el it u_previous

% Number of cells :
N_c = 6;

prjname_SOFC       = 'Solid Oxide Fuel Cell';
prjname_myProblem  = 'Real Time Optimization';

SOFC_data_nominal = data_SOFC_nominal(prjname_SOFC,N_c);
SOFC_data_plant = data_SOFC_plant(prjname_SOFC,N_c);
myProblem = data_myProblem(prjname_myProblem);

%%   Parameters
% Stream temperatures:
T_CH4_in  = 200.0; % [C]  methane stream  
T_H2O_in  = 200.0; % [C]  steam stream
T_air_in  = 30;    % [C]  air stream

T_in = [T_CH4_in T_H2O_in T_air_in] + SOFC_data_nominal.cst.K;

%% Initial guess for the optimization
% % Initial State :
T_r0    = 550; % [C]   fuel reformer

T_el0   = 650; % [C]   electrolyte
T_i0    = 650; % [C]   interconnect 
T_fuel0 = 650; % [C]   fuel channel
T_air0  = 650; % [C]   air channel 

T_b0    = 1000; % [C]  burner

T_h0    = 650; % [C]   heat exchanger - fuel side
T_c0    = 600; % [C]   heat exchanger - air side

T_0 = [T_r0 T_el0 T_i0 T_fuel0 T_air0 T_b0 T_h0 T_c0] + SOFC_data_nominal.cst.K;

% SCTP - CNTP
R = 8.314462175;    % [J/K/mol] 
P = 100e3;          % [Pa]
T_st = 273.15;      % [K]

q_CH4inp      = 0.376129940461854; % [L/min]     methane flow rate
q_AIRcathinp  = 17.2531441172150; % [L/min]     methane flow rate
Iinp          = 19;               % [A]         current

u_0 = [q_CH4inp,q_AIRcathinp,Iinp];

% Initial guess for the optimization problem
u0 = [u_0 T_0];

%% Constraints 
Pel_opt = [80 100 90];

% ub = [27.25,272.53/0.21, 30, 600*ones(1,1)+273.15 800*ones(1,4)+273.15 1600+273.15 1474 1474];
% lb = [1.36E-03,0.01/0.21, 0, 450*ones(1,6)+273.15   200+273.15 200+273.15];
ub = [27.25,272.53/0.21, 30, Inf*ones(1,8)];
lb = [1.36E-03,0.01/0.21, 0, -Inf*ones(1,8)];

Aeq = [];
beq = [];

Lair_upper = 10; 
Lair_lower = 3;
FU_upper   = 0.7;
%

kc = (6e+4)*N_c*R*T_st/(8*P*SOFC_data_nominal.cst.F);

A   = [2*Lair_lower,  -0.21,  0,  zeros(1,8);
      -2*Lair_upper,   0.21,  0,  zeros(1,8);
      -FU_upper,          0,  kc, zeros(1,8)];
B   =  [0;0;0];

%% -----------------------------------------------
% Optimal values of the plant at power setpoints
% ------------------------------------------------
myProblem.TC.TimeConstantOff = 1;
u_opt = [];
eff_opt = [];

for i = 1:3
    Ps_el = Pel_opt(i);
    [u_f,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,[],T_in,A,B,lb,ub,SOFC_data_plant,myProblem.OPT.options);
    [~,Ucell_opt(i),~,FU_opt(i),Lair_opt(i),~,~] = fPrimeMyProject(0,u_f(4:end),u_f(1:3),T_in,SOFC_data_nominal);
    u_opt = [u_opt; u_f];
    eff_opt = [eff_opt -FVAL];
end

%% ---------------------------------------
% Modifier adaptation optimization
% ----------------------------------------
% Optimal values
Pel_opt    = [80 90 100];
Ucell_opt  = 0.7;
eff_opt    = [0.4297   0.4260  0.4226];
inputs_opt = [0.3163   0.3589  0.402;
              12.4936 13.8194 15.1185;
              19.0476 21.4286 23.8095];

profile_setpoint = [2 2 2];  % profile setpoint

Ps_el = [];
Uopt_hist = [];
inputs_opt_hist = [];
eff_opt_hist = [];
for i = 1:size(profile_setpoint,2)
    Ps_el(i) = Pel_opt(profile_setpoint(i));
    inputs_opt_hist(:,i) = inputs_opt(:,profile_setpoint(i));
    eff_opt_hist(i) = eff_opt(profile_setpoint(i));
end

u_hist = []; % storing input values for all iterations

modif = zeros(3,4);
grad_toggle = 0; % activate gradient computation adaptation

Kca = [0 0.9 0.9];
Kgrad = 1;
deltaH = [1e-2 1e-2 1e-2];
u_previous = 0;
modif_ca_Ucell = [];
modif_ca_Pel   = [];
modif_gr_Ucell = [];
modif_gr_Pel   = [];
modif_gr_Effic = [];

for k = 1:size(Ps_el,2)
    it = k;
    [u_f,FVAL,EXITFLAG,OUTPUT,LAMBDA] = runobjconstr(u0,modif,T_in,A,B,lb,ub,SOFC_data_nominal,myProblem.OPT.options);
    
    [tempsteady_nominal] = OPTIM_SteadyState(u_f(1:3),u_f(4:end),T_in,SOFC_data_nominal);
    [~,Ucell_nomnl(k),Pel_nomnl(k),~,~,~,eta_sys_nomnl(k)] = fPrimeMyProject(0,tempsteady_nominal,u_f(1:3),T_in,SOFC_data_nominal);
    [tempsteady_plant] = OPTIM_SteadyState(u_f(1:3),u_f(4:end),T_in,SOFC_data_plant);
    [~,Ucell_plant(k),Pel_plant(k),~,~,~,eta_sys_plant(k)] = fPrimeMyProject(0,tempsteady_plant,u_f(1:3),T_in,SOFC_data_plant);
    
    modif(2,4) = (1-Kca(2))*modif(2,4) + Kca(2)*(Ucell_plant(k)-Ucell_nomnl(k));    % Cell voltage
    modif(3,4) = (1-Kca(3))*modif(3,4) + Kca(3)*(Pel_plant(k)-Pel_nomnl(k));        % Power set point
    
    if grad_toggle == 1
        for i = 1:3
%             u_nomsteady = u_f;
            u_grad   = u_f;
            u_grad(i) = u_grad(i) + deltaH(i);

            Tsteady_nomgrad    = OPTIM_SteadyState(u_grad,T_0,T_in,SOFC_data_nominal);
            u_grad(4:end)   = Tsteady_nomgrad;
            Tsteady_plagrad    = OPTIM_SteadyState(u_grad,T_0,T_in,SOFC_data_plant);
            u_grad(4:end)   = Tsteady_plagrad;
            
            [~,Ucell_nomnl1,Pel_nomnl1,~,~,~,eta_sys_nomnl1] = fPrimeMyProject(0,tempsteady_nominal,u_f(1:3),T_in,SOFC_data_nominal);
            [~,Ucell_nomnl2,Pel_nomnl2,~,~,~,eta_sys_nomnl2] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_nominal);
            [~,Ucell_plant1,Pel_plant1,~,~,~,eta_sys_plant1] = fPrimeMyProject(0,tempsteady_plant,u_f(1:3),T_in,SOFC_data_plant);
            [~,Ucell_plant2,Pel_plant2,~,~,~,eta_sys_plant2] = fPrimeMyProject(0,u_grad(4:end),u_grad(1:3),T_in,SOFC_data_plant);
            
            grad_Ucell_nomnl(i) = (Ucell_nomnl2-Ucell_nomnl1)/deltaH(i);
            grad_Ucell_plant(i) = (Ucell_plant2-Ucell_plant1)/deltaH(i);
            grad_Power_nomnl(i) = (Pel_nomnl2-Pel_nomnl1)/deltaH(i);
            grad_Power_plant(i) = (Pel_plant2-Pel_plant1)/deltaH(i);
            grad_Effic_nomnl(i) = (eta_sys_nomnl2-eta_sys_nomnl1)/deltaH(i);
            grad_Effic_plant(i) = (eta_sys_plant2-eta_sys_plant1)/deltaH(i);
        end
        modif(2,1:3) = (1-Kgrad)*modif(2,1:3) + Kgrad*(grad_Ucell_plant-grad_Ucell_nomnl); % Gradient cell voltage
        modif(3,1:3) = (1-Kgrad)*modif(3,1:3) + Kgrad*(grad_Power_plant-grad_Power_nomnl); % Gradient power
        modif(1,1:3) = (1-Kgrad)*modif(1,1:3) + Kgrad*(grad_Effic_plant-grad_Effic_nomnl); % Gradient efficiency
        
        modif_gr_Ucell = [modif_gr_Ucell; modif(2,1:3)];
        modif_gr_Pel   = [modif_gr_Pel;   modif(3,1:3)];
        modif_gr_Effic = [modif_gr_Effic; modif(1,1:3)];
    end
    
    modif_ca_Ucell = [modif_ca_Ucell; modif(2,4)];
    modif_ca_Pel   = [modif_ca_Pel; modif(3,4)];
    u0 = u_f;
    u_hist = [u_hist; u_f];
    u_previous = u_f(1:3);
end

%% ------------------------
% POST-PROCESSING
% -------------------------
% Optimal plant profile
% Pelopt = [Pel_opt(1)*ones(1,5), Pel_opt(2)*ones(1,5), Pel_opt(3)*ones(1,5)];
% etaopt = [eff_opt(1)*ones(1,5), eff_opt(2)*ones(1,5), eff_opt(3)*ones(1,5)];
% Ucellopt = [0.7*ones(1,5), 0.7*ones(1,5), 0.7*ones(1,5)];
% uopt = [u_opt(1,1:3)'*ones(1,5), u_opt(2,1:3)'*ones(1,5), u_opt(3,1:3)'*ones(1,5)];

% Optimal values
Popt_hist = Ps_el;
Uopt_hist = Ucell_opt*ones(1,size(Ps_el,2));

% Pelopt = [Pel_opt(1)*ones(1,5)]%, Pel_opt(2)*ones(1,5), Pel_opt(3)*ones(1,5)];
% etaopt = [eff_opt(1)*ones(1,5)]%, eff_opt(2)*ones(1,5), eff_opt(3)*ones(1,5)];
% Ucellopt = [0.7*ones(1,5)]%, 0.7*ones(1,5), 0.7*ones(1,5)];
% uopt = [u_opt(1,1:3)'*ones(1,5)]%, u_opt(2,1:3)'*ones(1,5), u_opt(3,1:3)'*ones(1,5)];
set(0,'DefaultFigureWindowStyle','docked')

figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([Pel_plant', Popt_hist']);
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex')
xlabel('iteration','Interpreter','latex')
ylabel('Power $P_{el}$ [W]','Interpreter','latex')
title('Power of the system','Interpreter','latex')
set(gca,'Box','off',...
'FontUnits','points',...
'FontWeight','normal','FontSize',12,...
'TickLabelInterpreter','latex')
grid on
%%
figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([eta_sys_plant', eff_opt_hist']);
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex')
xlabel('iteration','Interpreter','latex')
ylabel('efficiency $\eta$ [-]','Interpreter','latex')
title('Efficiency of the system','Interpreter','latex')
set(gca,'Box','off',...
'FontUnits','points',...
'FontWeight','normal','FontSize',12,...
'TickLabelInterpreter','latex')
grid on
%%
figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([Ucell_plant', Uopt_hist']);
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex','Location','southeast')
title('Voltage of the cell $U_{cell}$ [V]','Interpreter','latex')
xlabel('iteration','Interpreter','latex');
ylabel('$U_{cell}$ [V]','Interpreter','latex');
set(gca,'Box','off','FontUnits','points','FontWeight','normal','FontSize',11,...
'TickLabelInterpreter','latex')
grid on
%%
% figure
% set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
% stairs([modif_ca_Pel, modif_ca_Ucell]);
% legend({'power','voltage'},'Interpreter','latex')
% title('Constraint-values','Interpreter','latex')
% xlabel('iteration','Interpreter','latex');
% ylabel('Constraint-values [-]','Interpreter','latex');
% set(gca,'Box','off',...
% 'FontUnits','points',...
% 'FontWeight','normal','FontSize',11,...
% 'TickLabelInterpreter','latex')
% grid on
%%
figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([u_hist(:,1), inputs_opt_hist(1,:)']);
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex')
title('Flow rate of $\dot n_{CH_4}$ [L/min]','Interpreter','latex')
xlabel('iteration','Interpreter','latex');
ylabel('$\dot n_{CH_4}$ [L/min]','Interpreter','latex');
set(gca,'Box','off',...
'FontUnits','points',...
'FontWeight','normal','FontSize',11,...
'TickLabelInterpreter','latex')
grid on

%%
figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([u_hist(:,2), inputs_opt_hist(2,:)']);
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex')
title('Flow rate of air $\dot q_{Air}$ [L/min]','Interpreter','latex')
xlabel('iteration','Interpreter','latex');
ylabel('$\dot q_{Air}$ [L/min]','Interpreter','latex');
set(gca,'Box','off',...
'FontUnits','points',...
'FontWeight','normal','FontSize',11,...
'TickLabelInterpreter','latex')
grid on

%%
figure
set(gcf,'Color','w','Units','centimeters','Position',[30 10 16 9])
h = stairs([u_hist(:,3), inputs_opt_hist(3,:)']);
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(2).LineStyle = '--';
h(2).Color = 'k';
legend({'plant','optimal'},'Interpreter','latex')
title('Current of the fuel cell $I$ [A]','Interpreter','latex')
xlabel('iteration','Interpreter','latex');
ylabel('$I$ [A]','Interpreter','latex');
set(gca,'Box','off',...
'FontUnits','points',...
'FontWeight','normal','FontSize',11,...
'TickLabelInterpreter','latex')
grid on

%%
if isempty(modif_gr_Ucell) == 0
    figure
    stairs(modif_gr_Ucell)
    xlabel('iteration')
    ylabel('constraint-gradient modifier voltage')
    title('constraint-gradient modifier voltage')

    figure
    stairs(modif_gr_Pel)
    xlabel('iteration')
    ylabel('constraint-gradient modifier power')
    title('constraint-gradient modifier power')
    
    figure
    stairs(modif_gr_Effic)
    xlabel('iteration')
    ylabel('constraint-gradient modifier efficiency')
    title('constraint-gradient modifier efficiency')
end

%% ---------------------------
% Perturbation graphs
% ----------------------------
% tspan = [0:3:1000];
% [t_ode, y_ode] = ode15s(@(t, x) SIMUL_fPrimeMyProject(t,x,u_f,T_in,SOFC_data_plant), tspan, u_f(4:11), []);
% u_sim = u_f;
% for i = 1:size(y_ode,1)
%     u_sim = u_f;
%     u_sim(4:end) = y_ode(i,:);
%     
%     if t_ode(i)>200 && t_ode(i)<300
%         u_sim(1) = u_sim(1)+0.01;
%     elseif t_ode(i)>450 && t_ode(i)<550
%         u_sim(2) = u_sim(2)+1;
%     elseif t_ode(i)>700 && t_ode(i)<800
%         u_sim(3) = u_sim(3)+1;
%     else
%         u_sim(4:end) = y_ode(i,:);
%     end
%     [~,Ucell_plant_sim(i),Pel_plant_sim(i),~,~,~,eta_sys_plant_sim(i)] = SIMUL_fPrimeMyProject(t_ode(i),u_sim(4:end),u_sim(1:3),T_in,SOFC_data_plant);
% end
% figure
% stairs(t_ode,Pel_plant_sim)