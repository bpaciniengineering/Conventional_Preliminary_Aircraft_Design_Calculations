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
folder_name = uigetdir('C:\','Select Working Directory');
cd(folder_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial_Name  = 'Trial 1' ;
Description = ''        ;

M_cruise    = 0.85      ;
R           = 6500      ; %nm
AR          = 8.1       ; %assume about 8                       %ESTIMATE
e           = 0.8       ; %Oswald efficiency factor, assume 0.8 (Raymer 92)
tsfc        = 0.5       ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude_ci = 35000     ; %cruise altitude, ft
altitude_fi = 0000      ; %airfield alitude, ft
passengers  = 210       ; %persons
crew        = 5 + 3     ; %persons
baggage     = [40 4000] ; %lbs [allotment per person, additional cargo]
loiter_dur  = 0         ; %hrs

weight_max  = 1e6       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off

V_approach  = 150       ; %knots
V_stall = V_approach/1.3; %knots, based on approach speed estimate
Clmax_to    = 1.80      ; %assumed
Clmax_land  = 2.10      ; %assumed
L_takeoff   = 10500     ; %ft REQUIREMENT
L_landing= L_takeoff*0.6; %ft REQUIREMENT
M_climb     = M_cruise  ; %for now (see aircraft_mass.m)
rate_climb  = 2400      ; %ft/min
altitude_climbi = altitude_ci ; %ft, for now (see aircraft_mass.m)
theta_app   = 3.04      ; %approach angle, deg

% cruise parameters
C_D0_c      = 0.025      ; % assumed (at cruise)
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
    C_D0_c, C_DR_c, K1_c, K2_c, gamma, TR, g, carpet_x_lim, carpet_y_lim);


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
print('Carpet_Plot','-dpng')

%calculations
[ s_ref, thrust ] = TS_converter( wing_loading, thrust_loading, W_TO );

disp(sprintf('%0.0f Thrust (lbf)', thrust));
disp(sprintf('%0.0f Reference Area (ft^2)', s_ref));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TO TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fprintf(fid, sprintf('%0.0f Reference Area (ft^2)', s_ref));



