%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jan Bernhard & Bernardo Pacini                                        %%
%% MAE 332 - Aircraft Design                                             %%
%% Reserved Fuel Calculations                                            %%
%% Feb. 23, 2017 Thur                                                    %%
%% Modified: xx/xx/xxxx                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [weight_ratio] = fuel_reserve(AR,tsfc,Reserve_R,...
    M_cruise, altitude_ci, e, C_D0, C_DR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cruise speed
[airDens, airPres, temp, soundSpeed] = Atmos(altitude_ci);%kg/m^3 N/m^2 K m/s
soundSpeed = soundSpeed*2.23694;    %meter/sec -> mph
V_cruise = M_cruise*(soundSpeed);   %miles/hr
R = Reserve_R*1.15078;         %nm -> m
loiter_dur = 0.75;                  %hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Failed landing approach.
ratio_landing_1 = calculate_beta('land', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start through
ratio_startup = calculate_beta('takeoff', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise climb
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
ratio_landing_2 = calculate_beta('land', R, loiter_dur, 0, ...
    0, AR, e, C_D0+C_DR, tsfc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Ratio
weight_ratio = ratio_landing_1 * ratio_startup * ratio_climb * ratio_CO...
    * ratio_Loiter * ratio_landing_2;
end
