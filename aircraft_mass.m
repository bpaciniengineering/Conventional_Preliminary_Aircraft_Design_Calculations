%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 20, 2017 Mon                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_TO,W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, tsfc,...
    altitude, passengers, crew, baggage, loiter_dur, weight_max, graph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
altitude = altitude*0.3048;
[airDens, airPres, temp, soundSpeed] = Atmos(altitude);%kg/m^3 N/m^2 K m/s
% Convert values from SI to Imperial
airDens    = airDens * 0.0624;       %lb/ft^3
airPres    = airPres * 0.000145038;  %PSI
temp       = (9/5)*(temp - 273) + 32; %F
soundSpeed = soundSpeed*2.23694;     %convert to mph
% calculate cruise speed
V_cruise = M_cruise*(soundSpeed);
% calculating distance in miles
R = R*1.15078; %nm -> m
% Weight payload
W_pay = (passengers + crew)*(200) + (passengers + crew)*(baggage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Startup and Take-off
ratio_startup = 0.9725;                                         %ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration and Climb
ratio_climb = 1.0065 - 0.0325*M_cruise;                         %ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise Out
if M_cruise < 1
    L_D = AR + 10;                                              %M<1
else
    L_D = 11/sqrt(M);                                           %M>1
end
% Breguet Range Equation
% R = (V/tsfc) * (L_D) * ln(Wi/Wf) %lbfuel/h/lbt
ratio_CO = 1/(exp(R*((tsfc)/(V_cruise))/(L_D)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loiter
% E = (1/tsfc)*(L_D_max)*ln(Wi/Wf)
E = loiter_dur/(3600);
L_D_max = L_D / 0.94;
ratio_Loiter = 1/(exp((E * tsfc)/(L_D_max)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing
ratio_landing = 0.9725;                                         %ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Ratio
weight_ratio = ratio_startup* ratio_climb * ratio_CO * ...
    ratio_Loiter * ratio_landing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wto = linspace(1, weight_max, 1000);
w_fuel_mission = (wto - (wto*weight_ratio));
w_fuel = w_fuel_mission * 1.06;
raymer = wto.*(1-(1.02*wto.^(-0.06))) - w_fuel_mission - W_pay;
r = 0.*wto;
m = InterXmodified([wto;r],[wto;raymer]);
if graph == 1
    figure()
    hold on;
    plot(wto, wto.*(1-(1.02*wto.^(-0.06))) - w_fuel_mission- W_pay);
    plot(m(1,1), m(2,1),'r*');
    title('Raymer Equation')
    xlabel('Weight (lbs)');
    ylabel('N/A');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL WEIGHT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_TO = m(1,1);
W_fuel = (W_TO - (W_TO*weight_ratio));
W_empty = W_TO - W_fuel - W_pay;
end