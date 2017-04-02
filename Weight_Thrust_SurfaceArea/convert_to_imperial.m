%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nathan Wei                                                            %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% 7 March 2017                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function takes values from Atmos.m (kg/m3, Pa, K, m/s)
% and converts them into Imperial units (lbm/ft3, psi, degrees F, mph)
function [airDens_i, airPres_i, temp_i, soundSpeed_i] = ...
    convert_to_imperial(airDens, airPres, temp, soundSpeed)

airDens_i    = airDens * 0.0624;             %lbm/ft^3
airPres_i    = airPres / 6894.744;           %PSI
temp_i       = (9/5)*(temp - 273) + 32;      %F
soundSpeed_i = soundSpeed*2.23694;           %convert to mph

end

