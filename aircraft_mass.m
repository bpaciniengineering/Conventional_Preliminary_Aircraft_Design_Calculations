%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 23, 2017 Thur                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_TO,W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, tsfc,...
    altitude, passengers, crew, baggage, loiter_dur, weight_max, graph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
altitude = altitude*0.3048;
[airDens, airPres, temp, soundSpeed] = Atmos(altitude);%kg/m^3 N/m^2 K m/s
% Convert values from SI to Imperial
airDens    = airDens * 0.0624;          %lb/ft^3
airPres    = airPres * 0.000145038;     %PSI
temp       = (9/5)*(temp - 273) + 32;   %F
soundSpeed = soundSpeed*2.23694;        %convert to mph
% calculate cruise speed
V_cruise = M_cruise*(soundSpeed);
% calculating distance in miles
R = R*1.15078; %nm -> m
% Weight payload
if baggage(2) == 0
    W_pay = (passengers + crew)*(200) + (passengers + crew)*(baggage(1));
else
    W_pay = (passengers + crew)*(200) + (baggage(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Startup and Take-off
ratio_startup = 0.9725;                                         %ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration and Climb
if M_cruise < 1
    ratio_climb = 1.0065 - 0.0325*M_cruise;                        %ESTIMATE
else
    ratio_climb = 0.991 - 0.007*M_cruise - 0.01*M_cruise^2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise Out
if M_cruise < 1
    L_D = AR + 10;                                              %M<1
else
    L_D = 11/sqrt(M);                                           %M>1
end
% Breguet Range Equation
% R = (V/tsfc) * (L_D) * ln(Wi/Wf) %lbfuel/h/lbt
ratio_CO = 1/(exp(R*((tsfc)/(V_cruise))/(L_D)))
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
w = linspace(1, weight_max, 1000);
raymer = @(wto) wto*(1-(1.02*wto^(-0.06))) - ...
    (wto - (wto*weight_ratio))*1.06 - W_pay;
x0 = [1 weight_max];
W_TO = fzero(raymer, x0);
if graph == 1
    figure()
    hold on;
    plot(w, w.*(1-(1.02*w.^(-0.06))) - 1.06*(w-(w.*weight_ratio))-W_pay);
    plot(W_TO,0,'r*');
    title('Raymer Equation')
    xlabel('Weight (lbs)');
    ylabel('N/A');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL WEIGHT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_fuel = (W_TO - (W_TO*weight_ratio))/0.94; % 0.94 is trapped fuel
W_empty = W_TO - W_fuel - W_pay;
end
