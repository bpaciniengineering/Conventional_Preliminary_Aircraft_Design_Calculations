%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 27, 2017 Thur                                                    %%
%%                                                                       %%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file and a carpet plot with constraints        %%
%%                                                                       %%
%% Dependencies:| aircraft_mass.m | aircraft_surfacearea.m | Atmos.m |   %%
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

M_cruise    = 0.85      ;
R           = 2500      ; %nm
AR          = 8.1       ; %assume about 8                       %ESTIMATE
e           = 0.8       ; %Oswald efficiency factor, assume 0.8 (Raymer 92)
tsfc        = 0.5       ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude_ci = 35000     ; %cruise altitude, ft
altitude_fi = 0000      ; %airfield alitude, ft
passengers  = 8         ; %persons
crew        = 2         ; %persons
baggage     = [0 1000] ; %lbs [allotment per person, additional cargo]
loiter_dur  = 0         ; %hrs

weight_max  = 1e5       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off

V_approach  = 150       ; %knots
V_stall = V_approach/1.3; %knots, based on approach speed estimate
Clmax_to    = 1.80      ; %assumed
Clmax_land  = 2.10      ; %assumed
L_takeoff   = 10500     ; %ft REQUIREMENT
L_landing= L_takeoff*0.6; %ft REQUIREMENT
M_climb     = M_cruise*0.8; %for now (see aircraft_mass.m)
rate_climb  = 3500      ; %ft/min
altitude_climbi = 0     ; %ft, for now (see aircraft_mass.m)
theta_app   = 3.04      ; %approach angle, deg
n_max       = 6         ; % maximum load factor (FAR 23)
n_min       = -0.5      ; % minimum load factor (FAR 23)

% cruise parameters
C_D0_c      = 0.025     ; % assumed (at cruise)
C_DR_c      = 0         ; % assumed (clean configuration at cruise)
K1_c        = 1/(pi*AR*e); % induced drag correction factor
K2_c        = 0         ; % viscous drag correction factor
gamma       = 1.4       ; % specific heat ratio cp/cv, for air
TR          = 1         ; % assumed
g           = 32.174    ; %ft/s^2

carpet_x_lim = [50 150] ;
carpet_y_lim = [0 1]    ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Weight
[W_TO, W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, e, C_D0_c, ...
    C_DR_c, tsfc, altitude_ci, passengers, crew, baggage, loiter_dur,...
    weight_max, graph);

disp(sprintf('%0.0f Takeoff Weight (lbm)', W_TO)); 
disp(sprintf('%0.0f Fuel Weight (lbm)', W_fuel));
disp(sprintf('%0.0f Empty Weight (lbm)', W_empty));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
aircraft_surfacearea(M_cruise, R, AR, e, tsfc, altitude_ci, altitude_fi,...
    loiter_dur, altitude_climbi, V_approach, V_stall, Clmax_to,...
    Clmax_land, L_takeoff, L_landing, M_climb, rate_climb, theta_app,...
    C_D0_c, C_DR_c, K1_c, K2_c, gamma, TR, g, carpet_x_lim, carpet_y_lim, W_TO);


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
print(sprintf('Carpet Plot %s', Trial_Name),'-dpng')

%calculations
[ s_ref, thrust ] = TS_converter( wing_loading, thrust_loading, W_TO );

disp(sprintf('%0.0f Thrust (lbf)', thrust));
disp(sprintf('%0.0f Reference Area (ft^2)', s_ref));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(Trial_Name)
cd(Trial_Name)

text = sprintf('Preliminary Aircraft Design Calculations - %s.txt',...
    Trial_Name);
fid = fopen(text,'w');
fprintf(fid, sprintf('Preliminary Aircraft Design Calculations- %s \n',...
    Trial_Name));
fprintf(fid, sprintf('Description: %s \n \n',Description));

fprintf(fid, sprintf('%0.0f Takeoff Weight (lbm) \n', W_TO)); 
fprintf(fid, sprintf('%0.0f Fuel Weight (lbm) \n', W_fuel));
fprintf(fid, sprintf('%0.0f Empty Weight (lbm) \n \n', W_empty));
fprintf(fid, sprintf('%0.0f Thrust (lbf) \n', thrust));
fprintf(fid, sprintf('%0.0f Reference Area (ft^2) \n \n', s_ref));

fprintf(fid, sprintf('-------------------------------------------------'));
fprintf(fid, sprintf('\nInput Parameters: \n \n'));

fprintf(fid, sprintf('%0.2f Cruise Mach Number\n', M_cruise)); 
fprintf(fid, sprintf('%0.0f Range (Nm)\n', R)); 
fprintf(fid, sprintf('%0.1f Aspect Ratio\n', AR)); 
fprintf(fid, sprintf('%0.2f Oswald Efficiency Factor\n', e)); 
fprintf(fid, sprintf('%0.0f TSFC\n', tsfc)); 
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
