clear;close all; 

g = 9.81;       % sea level gravity acceleration [m/s²]
C_b = 0.47;     % form coefficent 
v =2 ;          % ship speed [m/s]
rho_H2O = 1000; % water density [kg/m³]
mu_H2O = 1e-3;  % water dynamic viscosity [Pa*s]

%% Power required for Drag

%%% Central Hull Data 

L_h = 25;        % length [m]
B_h = 3.46;      % width [m]
rho_ship = 2699; % ship material density [kg/m^3]

m = 45000;       % Ship mass [kg]
gamma = 10250;   % specific weight [N/m^3] 

V = m*g/gamma;   % draft volume [m³]

H_h = V /(L_h*B_h*C_b); % draft heigh [m]
S = 1.7*L_h*H_h+V/H_h;  % ship draft surface [m²]

% Estimation of the viscous drag introduced by the central hull

Re = rho_H2O*v*H_h/mu_H2O;     % Reynolds number 

K = 19*(V/(L_h^2*H_h))^2;
C_F = 0.075/(log10(Re)-2)^2;
C_V = C_F*(1+K);               % drag coefficient

R_v_c = 0.5*rho_H2O*v^2*C_V*S; % viscosity force [N]
P_h = R_v_c*v/0.7;             % Shaft Power central hull [W]
Results_central_hull = sprintf(' Viscous Drag Force central hull: %g kN\n Shaft Power central hull: %g MW',R_v_c*1e-3,P_h*1e-3);

%%% Side Hull Data

L_s = L_h/2;    % length [m]
H_s = 0.25;     % draft heigh [m]
B_s = 0.6 ;     % width [m]

V_s = L_s*H_s*B_s*C_b;        % draft volume (Archimede's principle) [m³]
S_s = 1.7*L_s*H_s+V_s/H_s;    % ship draft surface [m²]

% Estimation of the viscous drag introduced by the side hulls

Re_s = rho_H2O*v*H_s/mu_H2O;  % Reynolds number 

K_s = 19*(V_s/(L_s^2*H_s))^2;
C_F_s = 0.075/(log10(Re_s)-2)^2;
C_V_s = C_F_s*(1+K_s);        % drag coefficient

R_v_s = 0.5*rho_H2O*v^2*C_V_s*S_s;  % viscosity force [N]
P_s = R_v_s*v/0.7;                  % Shaft Power side hull [W]
Results_side_hull = sprintf('Viscous Drag Force side hull: %g kN\n Shaft Power side hull: %g MW',R_v_s*1e-3,P_s*1e-3);

R_v_tot = R_v_c+2*R_v_s;  % total viscosity force [N]
P_v_tot = P_h+2*P_s;      % total viscosity shaft Power [W]

Results_total_viscous = sprintf('Total Viscous Drag Force : %g kN\n Total Shaft Power: %g MW',R_v_tot*1e-3,P_v_tot*1e-3);

%%% Estimation of the wave drag 

C_w = 0.001;
R_w = 0.5*rho_H2O*C_w*v^2*(S+2*S_s);  % drag wave resistence 

R_tot = R_v_tot+R_w;                  % total drag force [N]
P_R = R_tot*v/0.7;                    % total shaft Power [W]
Results_total = sprintf('Total Drag Force : %g kN\nTotal Shaft Power for resistence: %g kW',R_tot*1e-3,P_R*1e-3);

%%% Estimation of net drag

Cn = 0.6; % Net Drag Coefficent
Sn = 0.9; % [m^2] Net twine area ( 4.5x0.4 m total net area)

R_n = 2*0.5*Cn*Sn*v^2*rho_H2O; % Drag produced by two nets

%% Total Power Required

P_tot = (P_R*1.25 + R_n*v);

%% Battery sizing

P_required = P_tot;         % Power required in the worst condition ( while catching )[W]
E_day = P_required*12;      % required energy from the batteries
DoD = 0.7;                  % Depth of discharge 
E_safe = E_day/DoD;         % Effective energy of the batteries
Volt = 48;                  % voltage of the selected battery
Amp = E_safe/Volt;          % capacity 
N = Amp/1000;               % Number of batteries

%% Solar Panel Power

P_density_sp = 150;         % power density of the selected solar panel
Sp_area = 192;              % Solar panel surface [m^2]
h_sun= 9;                   % Hours of sun 
P_p = Sp_area*P_density_sp; % Panels power [W]
E_p= P_p*h_sun;             % Ideal panels energy [Wh]
E_peff= E_p*0.8*0.93;       % Effective power energy [Wh]
E_rp= P_tot*12;             % Energy required for ship motion [Wh]
E_rptot = E_rp+E_day;       % Total energy required [Wh]

%% Wind Turbine Power

rho_aria = 1.226;              %[kg/m^3]
D = 1.6;                       % Turbine diameter [m]
A_in = pi*(D^2)/4;
C_be = 0.3;                    % Betz limit
v_w = 15;                      % Wind velocity [m/s]
P_e = 0.5*rho_aria*A_in*v_w^3; % Ideal turbine power [W]
P_dis = C_be*P_e*0.7;          % Effective turbine power [W]
E_dis = P_dis*24;              % Turbine energy [Wh]

%% Quantity of plastic 

q_min = 0.13;          % minimum density of plastic [kg/km]
q_max = 1.3;          % maximum density of plastic [kg/km]
velocity = v*3.6;     % velocity in km/h
plastic_day_min = q_min*velocity*24; 
plastic_day_max = q_max*velocity*24;
plastic_month_min = plastic_day_min*30;
plastic_month_max = plastic_day_max*30;
plastic_year_min = plastic_month_min*8;
plastic_year_max = plastic_month_max*8;

%% CO2 from Diesel engine
C_e = 692; %emission factor [g/kWh]
E_d = 167; % energy required for the motion during travel phase (day)[kWh]
E_n =73; % energy required for the motion during travel phase ( night)[kWh]
E_c = 200; % energy required for the motion during catchin phase [kWh]

g_CO2_departure_dayside = C_e*E_d; 
 
g_CO2_arrival_nightside = C_e*E_n; 

g_CO2_catching = C_e*E_c; 

g_CO2_catching_month = g_CO2_catching*30; % g of CO2 produced during one month of catching (30 days)
g_CO2_travel = (g_CO2_departure_dayside + g_CO2_arrival_nightside)*9; % g of CO2 produced during one travel (9 days)

%% structural mass estimation

s= 0.67;               % shape factor
k = 0.75;              % keel factor
V1 = s*k*L_h*H_h*B_h; % central hull full volume [m^3]
V_int = s*(L_h-0.06)*(H_h*0.75)*(B_h-0.06)-0.006*B_h*H_h*0.75*5; % central hull effective volume [m^3]
V2 = s*k*L_s*H_s*B_s;  % side hull volume [m^3]
B_sup= 4;              % beam of deadweight [m]
H_sup = 3;             % height of deadweight [m]
V_sup = B_sup*L_h*H_sup-(B_sup-0.06)*(L_h-0.06)*H_sup; % volume of deadweight
V_st = ((V1-V_int)+2*V2+V_sup)*1.5;                    % structural volume conservative estimation [m^3]

m_st = rho_ship*V_st;  % structural mass [kg]

%% internal systems mass

m_bat = 2000;          % mass of one battery [kg]
n= ceil(N);            % number of batteries
m_bat_tot = n*m_bat;   % total batteries mass [kg]

rho_p = 2.35;          % areal density solar panels [ kg/m^2]
m_p = rho_p*Sp_area;   % solar panels mass [kg]

m_m = 400;             % motor and net systems mass [kg]

%% total mass estimation

m_tot = m_st+m_m+m_p+m_bat_tot; % total mass [kg]
m_min = 40000;                  % minimum ship mass estimation [kg]
m_max = 60000;                  % maximum ship mass estimation [kg]

if m_min < m_tot < m_max
  disp( 'good mass estimtion')
 else
  disp('no good mass estimation')
end
