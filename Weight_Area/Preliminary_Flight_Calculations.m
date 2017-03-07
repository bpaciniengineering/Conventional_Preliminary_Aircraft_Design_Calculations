%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 27, 2017 Thur                                                    %%
%%                                                                       %%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file.                                          %%
%%                                                                       %%
%% Extra Dependencies: | aircraft_mass.m | atmos.m |                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SELECT WORKING DIRECTORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%folder_name = uigetdir('C:\','Select Working Directory');
%cd(folder_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial_Name  = 'Trial 1' ;
Description = ''        ;

M_cruise    = 0.8       ;
R           = 3500      ; %nm
AR          = 9         ; %assume about 8                       %ESTIMATE
e           = 0.8       ; %Oswald efficiency factor, assume 0.8 (Raymer 92)
tsfc        = 0.605     ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude_ci = 35000     ; %cruise altitude, ft
altitude_fi = 0000      ; % airfield alitude, ft
passengers  = 200       ; %persons
crew        = 10        ; %persons
baggage     = [50 0]    ; %lbs allotment passenger or crew
loiter_dur  = 1         ; %hrs

weight_max  = 1e6       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off

V_stall     = 137       ; %knots - better estimate possible?
V_approach  = 150       ; %knots
Clmax_to    = 1.80      ; %assumed
Clmax_land  = 2.10      ; %assumed
L_takeoff   = 10500     ; %ft REQUIREMENT
L_landing   = 3600      ; %ft REQUIREMENT
M_climb     = 0.7       ; %arbitrary
rate_climb  = 2000      ; %ft/min
altitude_climbi = 16000 ; %ft (guess)
theta_app   = 3         ; %approach angle, deg

% cruise parameters
C_D0_c      = 0.02      ; % assumed (at cruise)
C_DR_c      = 0         ; % assumed (clean configuration at cruise)
K1_c        = 1/(pi*AR*e); % induced drag correction factor
K2_c        = 0         ; % viscous drag correction factor
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
TR          = 1         ; % assumed
g           = 32.174    ; %ft/s^2

carpet_x_lim = [50 125];
carpet_y_lim = [0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Weight
[W_TO, W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, tsfc,...
    altitude_ci, passengers, crew, baggage, loiter_dur, weight_max, graph);

disp(sprintf('%0.0f Takeoff Weight', W_TO)); 
disp(sprintf('%0.0f Fuel Weight', W_fuel));
disp(sprintf('%0.0f Empty Weight', W_empty));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
%[] = aircraft_surfacearea();
altitude_c = altitude_ci*0.3048;
altitude_climb = altitude_climbi*0.3048;
altitude_f = altitude_fi*0.3048;
[airDens_c, airPres_c, temp_c, soundSpeed_c] = Atmos(altitude_c);%kg/m^3 N/m^2 K m/s
[airDens_f, airPres_f, temp_f, soundSpeed_f] = Atmos(altitude_f);
[airDens_climb, airPres_climb, temp_climb, soundSpeed_climb] = Atmos(altitude_climb);
[airDens_sl, airPres_sl, temp_sl, soundSpeed_sl] = Atmos(0);

% 84 degrees Fahrenheit field conditions (for midterm)
airDens_f       = 1.1644;   % kg/m^3
temp_f          = 28.889;   % degrees C
soundSpeed_f    = sqrt(1.4*287.1*(temp_f+273.15));

% Convert values from SI to Imperial
airDens_ci    = airDens_c * 0.0624;         %lbm/ft^3
airPres_ci    = airPres_c / 6894.744;    %PSI
temp_ci       = (9/5)*(temp_c - 273) + 32;  %F
soundSpeed_ci = soundSpeed_c*2.23694;       %convert to mph
airDens_fi    = airDens_f * 0.0624;         %lb/ft^3
airPres_fi    = airPres_f * 0.000145038;    %PSI
temp_fi       = (9/5)*(temp_f - 273) + 32;  %F
soundSpeed_fi = soundSpeed_f*2.23694;       %convert to mph
airDens_climbi    = airDens_climb * 0.0624;         %lb/ft^3
airPres_climbi    = airPres_climb * 0.000145038;    %PSI
temp_climbi       = (9/5)*(temp_climb - 273) + 32;  %F
soundSpeed_climbi = soundSpeed_climb*2.23694;       %convert to mph
airDens_sli = airDens_sl * 0.0624;          %air density at sea level
sigma = airDens_fi/airDens_sli;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stall
figure()
hax=axes; 
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
hold on;
V_stall = V_stall * 1.68781; %convert to ft/s
WS_stall = ((V_stall^2)*airDens_sli*Clmax_to)/(2*32.174);
% Plotted later for cosmetic reasons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off
WS = linspace(1,200);
TW_takeoff = ((20.9.*WS)/(sigma*Clmax_to)).*...
    (L_takeoff-69.6.*(WS./(sigma*Clmax_to)).^(.5)).^(-1);

plot(WS, TW_takeoff);
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 0 0]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Cruise Flight
% TW_CCF = (beta/alpha)*(k1*(beta/q)*WS + k2 + (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance
dHdt = rate_climb;
V_climb = (soundSpeed_climb*3.28084)*M_climb; % ft/s
beta_climb = 1.0065 - 0.0325*M_climb;
q_climb = 0.5*(airDens_climbi/g)*V_climb^2; % note: 1 lbm = 1/g slugs

% calculate alpha_tilde (for high bypass ratio turbofan engine) - MOVE
theta0 = (temp_climb/temp_sl)*(1+0.5*(gamma-1)*M_climb^2);
delta0 = (airPres_climb/airPres_sl)*...
    (1+0.5*(gamma-1)*M_climb^2)^(gamma/(gamma-1));

if theta0 <= TR
    alpha_climb = delta0*(1 - 0.49*M_cruise^0.5);
else
    alpha_climb = delta0*(1 - 0.49*M_cruise^0.5 - 3*(theta0-TR)/(1.5+M_cruise));
end

TW_climb = (beta_climb/alpha_climb).*(K1_c*(beta_climb/q_climb)*WS + K2_c + ...
    (C_D0_c + C_DR_c)./((beta_climb/q_climb)*WS) + (1/V_climb)*(dHdt));

%plot(WS, TW_climb, 'm'); % REQUIRED THRUST LOADING IS TOO HIGH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance

V_c = (soundSpeed_c*3.28084)*M_cruise; % ft/s
q = 0.5*(airDens_ci/g)*V_c^2; % note: 1 lbm = 1/g slugs

if M_cruise < 1
    L_D = AR + 10;
else 
    L_D = 11/sqrt(M_cruise);
end
% Breguet Range Equation
% R = (V/tsfc) * (L_D) * ln(Wi/Wf) %lbfuel/h/lbt
beta_c = 1/(exp(R*6076.12*((tsfc/3600)/(V_c))/(L_D)));

% calculate alpha_tilde (for high bypass ratio turbofan engine)
theta0 = (temp_c/temp_sl)*(1+0.5*(gamma-1)*M_cruise^2);
delta0 = (airPres_c/airPres_sl)*...
    (1+0.5*(gamma-1)*M_cruise^2)^(gamma/(gamma-1));

if theta0 <= TR
    alpha_c = delta0*(1 - 0.49*M_cruise^0.5);
else
    alpha_c = delta0*(1 - 0.49*M_cruise^0.5 - 3*(theta0-TR)/(1.5+M_cruise));
end

TW_cruise = (beta_c/alpha_c)*(K1_c*beta_c*WS/q + K2_c + ...
    (C_D0_c+C_DR_c)./(beta_c*WS/q));

plot(WS, TW_cruise, 'g');

%n = sqrt(1 + ((V^2))/g*R); % - here we assume n = 1 (no turning)

%TW_CLT = (beta/alpha)*(k1*n^2*(beta/q)*WS + k2*n + ...
  % (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loiter time in hours
%Estimate for (L/D)_max
LD_max = 0.5*sqrt(pi*e*AR/C_D0_c);

beta_l = 1/(exp((loiter_dur*tsfc)/LD_max));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing

WS_landing = (L_landing - 50/tan(deg2rad(theta_app)))*sigma*Clmax_land/79.4;
line([WS_landing WS_landing], get(hax,'YLim'),'Color',[0 1 1]);

%WS_landing = ((sigma * Clmax_land)/(79.4)) *(SL  - 50/tan(theta))

%plot(WS, TW_landing);
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 0 0]);
xlim(carpet_x_lim)
ylim(carpet_y_lim)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('Perliminary Aircraft Design Calculations - %s.txt',...
    Trial_Name);
fid = fopen(text,'w');
fprintf(fid, sprintf('Perliminary Aircraft Design Calculations- %s',...
    Trial_Name));
fprintf(fid, sprintf('\n')); 
fprintf(fid, sprintf('%s',Description));
fprintf(fid, sprintf('\n\n'));

fprintf(fid, sprintf('%0.0f Takeoff Weight', W_TO)); 
fprintf(fid, sprintf('\n')); 
fprintf(fid, sprintf('%0.0f Fuel Weight', W_fuel));
fprintf(fid, sprintf('\n')); 
fprintf(fid, sprintf('%0.0f Empty Weight', W_empty));
fprintf(fid, sprintf('\n')); 
