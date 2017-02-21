%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 20, 2017 Mon                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description: This code will output preliminary aircraft design        %%
%% calculations to a .txt file.                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dependencies: | aircraft_mass.m | atmos.m | InterXmodified.m |        %%
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
R           = 3600      ; %nm
AR          = 8         ; %assume about 8                       %ESTIMATE
tsfc        = 0.605       ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude    = 35000     ; %ft
passengers  = 200       ; %persons
crew        = 6         ; %persons
baggage     = 50        ; %lbs allotment passenger or crew
loiter_dur  = 0         ; %sec

V_stall     = 137       ; %knots
Clmax       = 1.5       ;      %assumed
L_runway    = 7000      ;    %assumed

weight_max  = 1e6       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY BELOW THIS POINT%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CALLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off Weight
[W_TO, W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, tsfc,...
    altitude, passengers, crew, baggage, loiter_dur, weight_max, graph);

disp(sprintf('%0.0f Takeoff Weight', W_TO)); 
disp(sprintf('%0.0f Fuel Weight', W_fuel));
disp(sprintf('%0.0f Empty Weight', W_empty));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Area
%[] = aircraft_surfacearea();
altitude = altitude*0.3048;
[airDens, airPres, temp, soundSpeed] = Atmos(altitude);%kg/m^3 N/m^2 K m/s
% Convert values from SI to Imperial
airDens    = airDens * 0.0624;       %lb/ft^3
airPres    = airPres * 0.000145038;  %PSI
temp       = (9/5)*(temp - 273) + 32; %F
soundSpeed = soundSpeed*2.23694;     %convert to mph

airDens_sl = 0.0765; %air density at sea level
sigma = airDens/airDens_sl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stall
figure()
hax=axes; 
hold on;
V_stall = V_stall * 1.68781; %convert to ft/s
WS_stall = ((V_stall^2)*airDens_sl*Clmax)/(2*32.174);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off
WS = linspace(1,200);
TW_takeoff = ((20.9.*WS)/(sigma*Clmax)).*...
    (L_runway-69.6.*(WS./(sigma*Clmax)).^(.5)).^(-1);



plot(WS, TW_takeoff);
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing



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
