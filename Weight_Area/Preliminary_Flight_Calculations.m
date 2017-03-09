%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 27, 2017 Thur                                                    %%
%%                                                                       %%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file.                                          %%
%%                                                                       %%
%% Extra Dependencies: | aircraft_mass.m | Atmos.m | calculate_alpha.m | %%
%%    calculate_beta.m | convert_to_imperial.m | dynamic_pressure.m      %%
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

M_cruise    = 0.85       ;
R           = 6500      ; %nm
AR          = 8         ; %assume about 8                       %ESTIMATE
e           = 0.8       ; %Oswald efficiency factor, assume 0.8 (Raymer 92)
tsfc        = 0.576     ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude_ci = 35000     ; %cruise altitude, ft
altitude_fi = 0         ; % airfield alitude, ft
passengers  = 200       ; %persons
crew        = 10        ; %persons
baggage     = [50 0]    ; %lbs allotment passenger or crew
loiter_dur  = 0         ; %hrs

weight_max  = 1e6       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off

V_approach  = 150       ; %knots
V_stall     = V_approach/1.3; %knots
Clmax_to    = 1.80      ; %assumed
Clmax_land  = 2.10      ; %assumed
L_takeoff   = 10500     ; %ft REQUIREMENT
L_landing   = 4800      ; %ft estimate
M_climb     = M_cruise  ; %for now (see aircraft_mass.m)
rate_climb  = 2000      ; %ft/min
altitude_climbi = altitude_ci ; %ft, for now (see aircraft_mass.m)
theta_app   = 3         ; %approach angle, deg

% cruise parameters
C_D0_c      = 0.02      ; % assumed (at cruise)
C_DR_c      = 0         ; % assumed (clean configuration at cruise)
K1_c        = 1/(pi*AR*e); % induced drag correction factor
K2_c        = 0         ; % viscous drag correction factor
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
TR          = 1         ; % assumed
g           = 32.174    ; %ft/s^2

carpet_x_lim = [50 150];
carpet_y_lim = [0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Weight
[W_TO, W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, e, C_D0_c, C_DR_c, tsfc,...
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

% convert nm to miles
R = R*1.15078;

% 84 degrees Fahrenheit field conditions (for midterm)
airDens_f       = 1.1644;   % kg/m^3
temp_f          = 28.889;   % degrees C
soundSpeed_f    = sqrt(1.4*287.1*(temp_f+273.15));

% Convert values from SI to Imperial
[airDens_ci, airPres_ci, temp_ci, soundSpeed_ci] = ...
    convert_to_imperial(airDens_c, airPres_c, temp_c, soundSpeed_c);
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = ...
    convert_to_imperial(airDens_f, airPres_f, temp_f, soundSpeed_f);
[airDens_climbi, airPres_climbi, temp_climbi, soundSpeed_climbi] = ...
    convert_to_imperial(airDens_climb, airPres_climb, temp_climb, soundSpeed_climb);
[airDens_sli, airPres_sli, temp_sli, soundSpeed_sli] = ...
    convert_to_imperial(airDens_sl, airPres_sl, temp_sl, soundSpeed_sl);

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
WS_stall = ((V_stall^2)*airDens_sli*Clmax_to)/(2*g);
% Plotted later for cosmetic reasons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off
WS = linspace(1,200);
TW_takeoff = ((20.9.*WS)/(sigma*Clmax_to)).*...
    (L_takeoff-69.6.*(WS./(sigma*Clmax_to)).^(.5)).^(-1);

plot(WS, TW_takeoff, 'b--');
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off - revised to include full take-off field length
index = (-28.43 + sqrt(28.43^2 - 4*(857.4-L_takeoff)*0.0185))/(2*0.0185);
% Take positive root only
% index2 = (-28.43 - sqrt(28.43^2 - 4*(857.4-L_takeoff)*0.0185))/(2*0.0185);

TW_takeoff1 = WS / (sigma*Clmax_to*index);
%TW_takeoff2 = WS / (sigma*Clmax_to*index2);
plot(WS, TW_takeoff1, 'b');
%plot(WS, TW_takeoff2, 'c--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Cruise Flight
% TW_CCF = (beta/alpha)*(k1*(beta/q)*WS + k2 + (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance
dHdt = rate_climb;
V_climb = (soundSpeed_climb*3.28084)*M_climb; % ft/s
beta_climb = calculate_beta('climb', R, loiter_dur, M_climb, ...
    soundSpeed_climbi*M_climb, AR, e, C_D0_c+C_DR_c, tsfc);
q_climb = dynamic_pressure(airDens_climbi, V_climb, g);

alpha_climb = calculate_alpha(temp_climb, airPres_climb, temp_sl, airPres_sl, ...
    gamma, M_climb, TR);

TW_climb = (beta_climb/alpha_climb).*(K1_c*(beta_climb/q_climb)*WS + K2_c + ...
    (C_D0_c + C_DR_c)./((beta_climb/q_climb)*WS) + (1/V_climb)*(dHdt));

%plot(WS, TW_climb, 'm'); % REQUIRED THRUST LOADING IS TOO HIGH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance

V_c = (soundSpeed_c*3.28084)*M_cruise; % ft/s
q_c = dynamic_pressure(airDens_ci, V_c, g);

beta_c = calculate_beta('cruise', R, loiter_dur, M_cruise, ...
    soundSpeed_ci*M_cruise, AR, e, C_D0_c+C_DR_c, tsfc);

alpha_c = calculate_alpha(temp_c, airPres_c, temp_sl, airPres_sl, ...
    gamma, M_cruise, TR);

TW_cruise = (beta_c/alpha_c)*(K1_c*beta_c*WS/q_c + K2_c + ...
    (C_D0_c+C_DR_c)./(beta_c*WS/q_c));

plot(WS, TW_cruise, 'g');

%n = sqrt(1 + ((V^2))/g*R); % - here we assume n = 1 (no turning)

%TW_CLT = (beta/alpha)*(k1*n^2*(beta/q)*WS + k2*n + ...
  % (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loiter time in hours
beta_l = calculate_beta('loiter', R, loiter_dur, 0, ...
    0, AR, e, C_D0_c+C_DR_c, tsfc);

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
text = sprintf('Preliminary Aircraft Design Calculations - %s.txt',...
    Trial_Name);
fid = fopen(text,'w');
fprintf(fid, sprintf('Preliminary Aircraft Design Calculations- %s',...
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
