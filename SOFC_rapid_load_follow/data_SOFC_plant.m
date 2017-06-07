function SOFC = data_SOFC_plant(projectName,N_c)
%% Data structure for the SOFC system variables and parameters

SOFC.name = projectName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Variables            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% State variables  -> x,dx
% input variables  -> u = [nCH4in, nO2in, nH2Oin, I]
% output variables -> y = [U_c P_el FU L_air epsilon_el epsilon_sys]


% States

% SOFC.xm  = zeros(1,4); % Fuel cell mass balance states
% SOFC.xt  = zeros(1,4); % Fuel cell temperature states
% SOFC.xr  = 0;          % Reformer temperature state
% SOFC.xh  = zeros(1,3); % Heat exchanger temperature states
% SOFC.xb  = 0;          % Burner temperature state

% State derivatives

% SOFC.dxm  = zeros(1,4); 
% SOFC.dxt  = zeros(1,4);
% SOFC.dxr  = 0;          
% SOFC.dxh  = zeros(1,3); 
% SOFC.dxb  = 0; 

% Inputs and outputs

% SOFC.u  = [];
% SOFC.y  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Fuel cell parameters      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Anode         -> AN
%  Cathode       -> CA
%  Electrolyte   -> EL
%  Interconector -> IC
%  Channel       -> CH
%  Cell          -> CE

% Anode parameters

SOFC.AN.h   = 5e-4;  % [m]
SOFC.AN.l_x = 0.15;  % [m]
SOFC.AN.l_y = 0.072; % [m]
SOFC.AN.rho = 7.5*6600;  % [kg/m3]
SOFC.AN.C_p = 400.0; % [J/kg/K]
SOFC.AN.e_a = 0.9;   % [-]

% Cathode parameters

SOFC.CA.h   = 5e-5;  % [m]
SOFC.CA.l_x = 0.15;  % [m]
SOFC.CA.l_y = 0.072; % [m]
SOFC.CA.rho = 7.5*6600;  % [kg/m3]
SOFC.CA.C_p = 400.0; % [J/kg/K]
SOFC.CA.e_c = 0.9;   % [-]

% Electrode parameters

SOFC.EL.h   = 1e-6;  % [m]
SOFC.EL.l_x = 0.15;  % [m]
SOFC.EL.l_y = 0.072; % [m]
SOFC.EL.rho = 7.5*6600;  % [kg/m3]
SOFC.EL.C_p = 400.0; % [J/kg/K]
SOFC.EL.e_e = 0.9;   % [-]

% Interconnect parameters

SOFC.IC.h   = 0.001; % [m]
SOFC.IC.l_x = 0.15;  % [m]
SOFC.IC.l_y = 0.072; % [m]
%SOFC.IC.rho = 6110;  % [kg/m3]
SOFC.IC.C_p = 400.0; % [J/kg/K]
SOFC.IC.e_i = 0.9;   % [-]
SOFC.IC.R_MIC = 0.00987214; % [Ohm]
SOFC.IC.rho = 7.5*6110;  % [kg/m3]
% Channel parameters

SOFC.CH.h_a  = 0.00075;     %[m]
SOFC.CH.h_c  = 0.001;       %[m]
SOFC.CH.w_a  = 0.0035;      %[m]
SOFC.CH.w_c  = 0.0035;      %[m]

SOFC.CH.A_a  = 0.0025;      %[m2] ?? estimation
SOFC.CH.A_c  = 0.0025;      %[m2] ?? estimation
SOFC.CH.V_a  = 0.2*5/384;         %[m3] ?? estimation
SOFC.CH.V_c  = 0.2*5/384;         %[m3] ?? estimation
SOFC.CH.C_fa = 0.7500;      %[-]  ?? estimation
SOFC.CH.C_fc = 0.7500;      %[-]  ?? estimation

% SOFC.CH.T_a  = zeros(1,2);  % [K] , T_a(1) = T_a_in
% SOFC.CH.T_c  = zeros(1,2);  % [K] , T_c(2) = T_c_out
% SOFC.CH.T_e  = 0;           % [K]
% SOFC.CH.T_i  = 0;           % [K]

% SOFC.CH.nH2  = zeros(1,3);  % [mol/s]
% SOFC.CH.nH2O = zeros(1,3);  % [mol/s]
% SOFC.CH.nO2  = zeros(1,3);  % [mol/s]
% SOFC.CH.nN2  = zeros(1,3);  % [mol/s]
% SOFC.CH.nCO2 = zeros(1,3);  % [mol/s]
% SOFC.CH.nCH4 = zeros(1,3);  % [mol/s]
% SOFC.CH.nCO  = zeros(1,3);  % [mol/s]

SOFC.CH.r_H2O = 0.03092783478941;
SOFC.CH.r_N2  = 3.76190444967126;

% Cell parameters

SOFC.CE.N_c    = N_c;    % [-]
SOFC.CE.h_fl   = 0.01;   % [m]
SOFC.CE.l_xfl  = 0.192;  % [m]
SOFC.CE.l_yfl  = 0.072;  % [m]
SOFC.CE.l_insl = 0.03;   % [m]
SOFC.CE.A_loss = 0.0;    % [m2]
SOFC.CE.A_act  = 60e-4;  % [m2]
SOFC.CE.C_p    = 740.0;  % [J/kg/K]   ?? estimation
SOFC.CE.rho    = 4600.0; % [kg/Cell]  ?? estimation
SOFC.CE.FU_adj = 0.15;

SOFC.CE.e_st   = 0.8; % [-]
SOFC.CE.e_in   = 0.8; % [-]
SOFC.CE.e_VF   = 2/3; % [-]
% SOFC.CE.T_env  = 0.0; % [K]
% SOFC.CE.T_f    = 0.0; % [K]

% Kinetic parameters

SOFC.kinetic.E_actcath  = 153260.5;
SOFC.kinetic.k0         = 4.103096e11;
SOFC.kinetic.E_disscath = 1.543785;
SOFC.kinetic.R0         = 9.225228e-14; 
SOFC.kinetic.E_el       = 7.622e4;
SOFC.kinetic.sig0_el    = 16300;


% SOFC.kinetic.E_actcath  = 150000.0;
% SOFC.kinetic.k0         = 4.5e11;
% SOFC.kinetic.E_disscath = 1.54; 
% SOFC.kinetic.R0         = 10e-14; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Reformer parameters       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions propagation

% SOFC.ref.r0            = 0.001; 
% SOFC.ref.ksi0          = 0.0005;
% SOFC.ref.nCH4in        = 0;

% Physical properties

SOFC.ref.rho           = 1000;     % [kg/m^3]
SOFC.ref.Cp            = 1000;     % [J/kg/K]
SOFC.ref.Vr            = 0.00012;  % [m^3]
SOFC.ref.R             = 0.006;    % [m]

% Temperatures

% SOFC.ref.TCH4in        = 0;        % [K]
% SOFC.ref.TH2Oin        = 0;        % [K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Heat exchanger parameters    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOFC.hex.Tc_in         = 0;       % [K]

% Dimensions

SOFC.hex.L             = 0.3;     % [m]
SOFC.hex.W             = 0.15;    % [m]
SOFC.hex.H             = 0.005;    % [m]
SOFC.hex.e             = 0.0005;   % [m]
SOFC.hex.nPlate        = 10;      % [m]

% Metal characteristics

SOFC.hex.CpMet        = 0.49e3;   % [J/K/kg]
SOFC.hex.rhoMet       = 7750;     % [kg/m^3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Burner parameters         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOFC.burn.V    = 0.1;            % [m^3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Physical constants        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOFC.cst.p_atm = 101325.0000;    % [Pa]
SOFC.cst.K     = 273.15;         % [K]
SOFC.cst.R     = 8.314462175;    % [J/K/mol]     
SOFC.cst.F     = 96486.00000;    % [C/mol]
SOFC.cst.kB    = 1.3806488e-23;  % [J/K]
SOFC.cst.sigma = 5.670373e-8;    % [W/m2/k4]
SOFC.cst.T_ref  = 298.15;        % [K]

% Molar masses

SOFC.cst.MH2    =  2.016e-3;      % [kg/mol]
SOFC.cst.MH2O   = 18.016e-3;      % [kg/mol]
SOFC.cst.MO2    = 31.998e-3;      % [kg/mol]
SOFC.cst.MN2    = 28.014e-3;      % [kg/mol]   
SOFC.cst.MCH4   = 16.040e-3;      % [kg/mol]
SOFC.cst.MCO    = 28.0101e-3;     % [kg/mol]   
SOFC.cst.MCO2   = 44.0095e-3;     % [kg/mol]

% Standard enthalpy

SOFC.cst.H_H2    = 0;             % [J/mol]
SOFC.cst.H_H2O   = -241.83e3;     % [J/mol]
SOFC.cst.H_O2    = 0;             % [J/mol]
SOFC.cst.H_N2    = 0;             % [J/mol]   
SOFC.cst.H_CO    = -110.53e3;     % [J/mol]
SOFC.cst.H_CO2   = -393.52e3;     % [J/mol]
SOFC.cst.H_CH4   = -73.4e3;       % [J/mol]

% LHV

SOFC.cst.LHVCH4 = 50.009e+6;      % [J/kg]
SOFC.cst.LHV    = 119.96e+6;      % [J/kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% UPDATE STRUCTURE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c = SOFC.CE.N_c;

h_an   = SOFC.AN.h;
h_ca   = SOFC.CA.h;
h_el   = SOFC.EL.h;
h_MIC  = SOFC.IC.h;
h_fuel = SOFC.CH.h_a;
h_air  = SOFC.CH.h_c;
h_fl   = SOFC.CE.h_fl;
l_xfl  = SOFC.CE.l_xfl;
l_yfl  = SOFC.CE.l_yfl;

h_RE = h_an + h_ca + h_el + h_MIC + h_fuel + h_air;

A_TB    = l_xfl*l_yfl;
A_front = (h_RE*N_c + 2*h_fl)*l_xfl;
A_side = (h_RE*N_c + 2*h_fl)*l_yfl; 

SOFC.CE.A_loss = 2*(A_TB + A_front + A_side);

end

