%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 23, 2017 Thur                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_TO,W_fuel, W_empty] = aircraft_mass(M_cruise, R, AR, e, C_D0, C_DR, tsfc,...
    altitude, passengers, crew, baggage, loiter_dur, weight_max, graph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
altitude = altitude*0.3048;
[airDens, airPres, temp, soundSpeed] = Atmos(altitude);%kg/m^3 N/m^2 K m/s
% Convert values from SI to Imperial
[airDens, airPres, temp, soundSpeed] = ...
    convert_to_imperial(airDens, airPres, temp, soundSpeed);
% calculate cruise speed
V_cruise = M_cruise*(soundSpeed);
%convert from nm -> m
R = R*1.15078;
% convert from seconds to hours
E = loiter_dur/(3600);


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
ratio_startup = calculate_beta('takeoff', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration and Climb
ratio_climb = calculate_beta('climb', R, loiter_dur, M_cruise, ...
    V_cruise, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise Out
ratio_CO = calculate_beta('cruise', R, loiter_dur, M_cruise, ...
    V_cruise, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loiter
ratio_Loiter = calculate_beta('loiter', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing
ratio_landing = calculate_beta('land', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Ratio
weight_ratio = ratio_startup* ratio_climb * ratio_CO * ...
    ratio_Loiter * ratio_landing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = linspace(1, weight_max, 1000);
raymer = @(wto) wto*(1-(1.02*wto^(-0.06))) - ...
    (wto - (wto*weight_ratio))/0.94 - W_pay;
x0 = [1 weight_max];
W_TO = fzero(raymer, x0);
if graph == 1
    figure()
    hold on;
    plot(w, w.*(1-(1.02*w.^(-0.06))) - (w-(w.*weight_ratio))/0.94-W_pay);
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
