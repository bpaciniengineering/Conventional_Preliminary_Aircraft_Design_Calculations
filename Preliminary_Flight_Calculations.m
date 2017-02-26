%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 23, 2017 Thur                                                    %%
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
folder_name = uigetdir('C:\','Select Working Directory');
cd(folder_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial_Name  = 'Trial 1' ;
Description = ''        ;

M_cruise    = 0.85      ; 
R           = 2500      ; %nm
AR          = 8         ; %assume about 8                       %ESTIMATE
tsfc        = 0.7       ; %0.45<=tsfc<=1.2 - check engine manufacturer
altitude    = 35000     ; %ft
passengers  = 210       ; %persons
crew        = 0         ; %persons
baggage     = [4000 1]  ; %lbs, 0-crew/passenger allotment/1-total payload
loiter_dur  = 0         ; %sec

weight_max  = 1e6       ; %max of weight range
graph       = 1         ; %1/0 for plot on/off

V_stall     = 137       ; %knots
V_approach  = 150       ; %knots
Clmax       = 1.5       ; %assumed
L_takeoff   = 10500     ; %ft REQUIREMENT
L_landing   = 3600      ; %ft REQUIREMENT
rate_climb  = 3500      ; %ft/min
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alt = altitude*0.3048;
[airDens, airPres, temp, soundSpeed] = Atmos(alt);%kg/m^3 N/m^2 K m/s
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
    (L_takeoff-69.6.*(WS./(sigma*Clmax)).^(.5)).^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Cruise Flight
% beta = 
% alpha =
% k1 = 
% k2 = 
q = dynamic_viscosity(alt);

% TW_CCF = (beta/alpha)*(k1*(beta/q)*WS + k2 + (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance
% beta = 
% alpha =
% k1 = 
% k2 = 
q = dynamic_viscosity(alt);
dHdt = rate_climb;

%TW_CP = (beta/alpha)*(k1*(beta/q)*WS + k2 + (CD_O + CD_R)/((beta/q)*WS)...
 %   + (1/V)*(dHdt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance
% beta = 
% alpha =
% k1 = 
% k2 = 
% R = 
g = 32; %ft/s^2
q = dynamic_viscosity(alt);
%n = sqrt(1 + ((V^2))/g*R);

%TW_CLT = (beta/alpha)*(k1*n^2*(beta/q)*WS + k2*n + ...
  % (CD_O + CD_R)/((beta/q)*WS));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing
%SL = 
%theta = 

%WS_landing = ((sigma * Clmax)/(79.4)) *(SL  - 50/tan(theta))

plot(WS, TW_takeoff);
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 0 0]);


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
