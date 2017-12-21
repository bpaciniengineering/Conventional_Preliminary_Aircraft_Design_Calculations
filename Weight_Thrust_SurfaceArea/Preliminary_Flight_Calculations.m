%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini and Nathan Wei and Raj Balaji and Leif Fredericks     %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 27, 2017 Thur                                                    %%
%%                                                                       %%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file and a carpet plot with constraints        %%
%%                                                                       %%
%% Dependencies:| aircraft_weight.m | aircraft_carpetplot.m | Atmos.m |  %%
%%              | calculate_alpha.m | calculate_beta.m | TS_converter.m |%% 
%%              | convert_to_imperial.m | dynamic_pressure.m |           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial_Name  = 'Trial 1' ;
Description = ''        ;
% ON/OFF SELECTORS
output      = 1         ; %1/0 for output files
weight      = 1         ; %1/0 for weight output, generally set to 1
carpet_plot = 1         ; %1/0 for carpet plot
TO_landing_lengths = 1  ; %1/0 - requires weight & carpet_plot 
box_wing = 1            ; % 1/0 [box wing aircraft]/[conventional aircraft]
graph       = 0         ; %1/0 for graph of weight equation solution

% PERFORMANCE ESTIMATES
AR          = 8.1       ; % Aspect ratio: ESTIMATE then ITERATE
LDC         = 15        ; % L/D at cruise conditions (based on design C_L),
                          % ESTIMATE (AR + 10 if M<1; 11/sqrt(M)
                          % if M>1) then ITERATE
LD          = 18        ; % L/D max: ESTIMATE (LDC / 0.94) then ITERATE
h_b         = 0.22      ; % height to span ratio (height is vertical 
                          % distance between wings) for BOX WING, else 0
e           = 0.8       ; % Oswald efficiency factor, assume 0.8(Raymer 92)

% do not modify this conditional
if box_wing == 1
    e   = 0.8 * (0.44 + 2.219* (h_b)) / (0.44 + 0.9594*(h_b));
             % http://www.fzt.haw-hamburg.de/pers/Scholz/Airport2030/
             % Airport2030_PUB_DLRK_11-09-27.pdf
end

tsfc        = 0.5       ; % 0.45<=tsfc<=1.2 - check engine manufacturer
Clmax_to    = 1.80      ; % assumed: ESTIMATE then ITERATE
Clmax_land  = 2.10      ; % assumed: ESTIMATE then ITERATE
C_D0_c      = 0.025     ; % assumed (at cruise): ESTIMATE then ITERATE
C_DR_c      = 0         ; % assumed (clean configuration at cruise)
K1_c        = 1/(pi*AR*e); % induced drag correction factor, 
   % see Raymer chapter 12 for more rigorous estimation methods, esp M>0.5 
K2_c        = 0         ; % viscous drag correction factor

ws          = 75        ; % wing loading (to check design point Landing FL)

% AIRCRAFT REQUIREMENTS 
M_cruise    = 0.89      ;
R           = 2500      ; %nm
Reserve_R   = 100       ; %nm
altitude_ci = 35000     ; %cruise altitude, ft
altitude_fi = 0000      ; %airfield alitude, ft
passengers  = 8         ; %persons
crew        = 2         ; %persons
baggage     = [0 1000]  ; %lbs [allotment per person, additional cargo]
loiter_dur  = 0         ; %hrs

weight_max  = 1e5       ; %max of weight range order of magnitude

% FLIGHT REQUIREMENTS
V_approach  = 150       ; %knots
V_stall = V_approach/1.3; %knots, based on approach speed estimate
L_takeoff   = 4000      ; %ft REQUIREMENT
L_landing= 3600         ; %ft REQUIREMENT
M_climb     = M_cruise*0.8; %for now (see aircraft_mass.m)
rate_climb  = 3500      ; %ft/min
altitude_climbi = 0     ; %ft, for now (see aircraft_mass.m)
theta_app   = 3.04      ; %approach angle, deg
n_max       = 6         ; % maximum load factor (FAR 23)
n_min       = -0.5      ; % minimum load factor (FAR 23)

% CRUISE PARAMETERS
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
TR          = 1         ; % assumed
g           = 32.174    ; %ft/s^2

% CARPET PLOT LIMITS
carpet_x_lim = [50 100] ; % W/S ratio, generally between 50 and 100
carpet_y_lim = [0 1]    ; % T/W ratio, generally between .4 and .7

% ENGINE QUALITIES
Thrust      = 10000;   %lbf/engine (to check design point TOFL)
Number_Engines = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LANDING FIELD LENGTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
altitude_c = altitude_ci*0.3048;
altitude_climb = altitude_climbi*0.3048;
altitude_f = altitude_fi*0.3048;
[airDens_c, airPres_c, temp_c, soundSpeed_c] = Atmos(altitude_c); % SI
[airDens_f, airPres_f, temp_f, soundSpeed_f] = Atmos(altitude_f);
[airDens_climb, airPres_climb, temp_climb, soundSpeed_climb] = ...
    Atmos(altitude_climb);
[airDens_sl, airPres_sl, temp_sl, soundSpeed_sl] = Atmos(0);

% Convert values from SI to Imperial
[airDens_ci, airPres_ci, temp_ci, soundSpeed_ci] = ...
    convert_to_imperial(airDens_c, airPres_c, temp_c, soundSpeed_c);
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = ...
    convert_to_imperial(airDens_f, airPres_f, temp_f, soundSpeed_f);
[airDens_climbi, airPres_climbi, temp_climbi, soundSpeed_climbi] = ...
    convert_to_imperial(airDens_climb, airPres_climb, temp_climb, ...
    soundSpeed_climb);
[airDens_sli, airPres_sli, temp_sli, soundSpeed_sli] = ...
    convert_to_imperial(airDens_sl, airPres_sl, temp_sl, soundSpeed_sl);

sigma = airDens_fi/airDens_sli;

Landing_Dist = 79.4 * ws/(sigma * Clmax_land) + 50/tan(deg2rad(theta_app));
disp(sprintf('%0.0f Landing Field Length (ft)', Landing_Dist)); 

%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Weight
if weight == 1
[W_TO, W_fuel, W_empty] = aircraft_weight(M_cruise, R, AR, e, C_D0_c, ...
    C_DR_c, tsfc, altitude_ci, passengers, crew, baggage, loiter_dur,...
    Reserve_R, weight_max, graph, LD, LDC);

disp(sprintf('%0.0f Takeoff Weight (lbm)', W_TO)); 
disp(sprintf('%0.0f Fuel Weight (lbm)', W_fuel));
disp(sprintf('%0.0f Empty Weight (lbm)', W_empty));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
if carpet_plot == 1
try
aircraft_carpetplot(M_cruise, R, AR, e, tsfc, altitude_ci, altitude_fi,...
    loiter_dur, altitude_climbi, V_approach, V_stall, Clmax_to,...
    Clmax_land, L_takeoff, L_landing, M_climb, rate_climb, theta_app,...
    C_D0_c, C_DR_c, K1_c, K2_c, gamma, TR, g, carpet_x_lim, ...
    carpet_y_lim, W_TO, LD, LDC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THRUST AND SURFACE AREA VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Thrust Loading [T_0/W_g]:','Wing Loading [W_g/S]:'};
dlg_title = 'User Input';
num_lines = 1;
defaultans = {'Value...','Value...'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

thrust_loading = str2num(answer{1});
wing_loading = str2num(answer{2});
plot(wing_loading,thrust_loading,'b.', 'MarkerSize', 20);

%calculations
[ s_ref, thrust ] = TS_converter( wing_loading, thrust_loading, W_TO );

disp(sprintf('%0.0f Thrust (lbf)', thrust));
disp(sprintf('%0.0f Reference Area (ft^2)', s_ref));
catch
    errordlg('Please enable aircraft weight');
    error('Please enable aircraft weight');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAKEOFF DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
if TO_landing_lengths == 1
[ TOFL , index] = TO_distance(V_stall, Clmax_to, s_ref, ... 
    Thrust, W_TO, Number_Engines, altitude_fi);
disp(sprintf('%0.0f Required Takeoff Eield Length (ft)', TOFL));
end
catch
    errordlg('Please enable both weight and Carpet Plot calculations')
    error('Please enable both weight and Carpet Plot calculations')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output == 1
folder_name = uigetdir('C:\','Select Working Directory');
mkdir(sprintf('%s/%s', folder_name, Trial_Name));

text = sprintf('Preliminary Aircraft Design Calculations - %s.txt',...
    Trial_Name);
fid = fopen(sprintf('%s/%s/%s',folder_name, Trial_Name, text),'w');
fprintf(fid, sprintf('Preliminary Aircraft Design Calculations- %s \n',...
    Trial_Name));
fprintf(fid, sprintf('Description: %s \n \n',Description));
if weight == 1
fprintf(fid, sprintf('%0.0f Takeoff Weight (lbm) \n', W_TO)); 
fprintf(fid, sprintf('%0.0f Fuel Weight (lbm) \n', W_fuel));
fprintf(fid, sprintf('%0.0f Empty Weight (lbm) \n \n', W_empty));
end
if carpet_plot == 1
fprintf(fid, sprintf('%0.0f Thrust (lbf) \n', thrust));
fprintf(fid, sprintf('%0.0f Reference Area (ft^2) \n \n', s_ref));
end
if TO_landing_lengths == 1
fprintf(fid, sprintf('%0.0f Takeoff Field Length (ft) \n \n',TOFL));
fprintf(fid,sprintf('%0.0f Landing Field Length (ft) \n \n',Landing_Dist));
end

fprintf(fid, sprintf('------------------------------------------------'));
fprintf(fid, sprintf('\nInput Parameters: \n \n'));

fprintf(fid, sprintf('%0.2f Cruise Mach Number\n', M_cruise)); 
fprintf(fid, sprintf('%0.0f Range (Nm)\n', R)); 
fprintf(fid, sprintf('%0.1f Aspect Ratio\n', AR)); 
fprintf(fid, sprintf('%0.2f Oswald Efficiency Factor\n', e)); 
fprintf(fid, sprintf('%0.3f TSFC\n', tsfc)); 
fprintf(fid, sprintf('%0.0f Altitude_ci (ft) \n', altitude_ci)); 
fprintf(fid, sprintf('%0.0f Altitude_fi (ft) \n', altitude_fi)); 
fprintf(fid, sprintf('%0.0f Passengers \n', passengers)); 
fprintf(fid, sprintf('%0.0f Crew \n', crew)); 
fprintf(fid, sprintf('%0.0f Baggage (lbm) \n', baggage)); 
fprintf(fid, sprintf('%0.2f Loiter Duration (hrs) \n \n', loiter_dur)); 

fprintf(fid, sprintf('%0.1f Approach Velocity (knots) \n', V_approach)); 
fprintf(fid, sprintf('%0.2f Stall Velocity (knots) \n', V_stall)); 
fprintf(fid, sprintf('%0.2f Takeoff CL_max \n', Clmax_to)); 
fprintf(fid, sprintf('%0.2f Land CL_max \n', Clmax_land)); 
fprintf(fid, sprintf('%0.2f Climb Mach Number \n', M_climb)); 
fprintf(fid, sprintf('%0.2f Rate Climb (ft/min) \n', rate_climb)); 
fprintf(fid, sprintf('%0.2f Climb Altitude (ft) \n', altitude_climbi)); 
fprintf(fid, sprintf('%0.2f Approach angle [theta] (deg)\n\n', theta_app)); 

fprintf(fid, sprintf('%0.3f Coefficient of Drag \n', C_D0_c)); 
fprintf(fid, sprintf('%0.3f Coefficient of Residual Drag \n', C_DR_c)); 
fprintf(fid, sprintf('%0.3f Induced Drag Correction [K1] \n', K1_c)); 
fprintf(fid, sprintf('%0.3f Viscous Drag Correction [K2] \n', K2_c)); 
fprintf(fid, sprintf('%0.3f Specific Heat Ratio [gamma] \n', gamma)); 
fprintf(fid, sprintf('%0.1f TR \n', TR)); 
fprintf(fid, sprintf('%0.1f Gravitational Constant (ft/s^2) \n', g)); 

if carpet_plot == 1
print(sprintf('%s/%s/Carpet Plot %s', folder_name, Trial_Name, ...
Trial_Name),'-dpng')
end
end