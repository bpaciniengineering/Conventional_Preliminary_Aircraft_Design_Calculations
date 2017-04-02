%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jan Bernhard & Bernardo Pacini                                        %%
%% MAE 332 - Aircraft Design                                             %%
%% Reserved Fuel Calculations                                            %%
%% Mar. 08, 2017                                                         %%
%% Modified: xx/xx/xxxx                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE:  Estimates the take off field length with regards to number of
%       engines 

function [ TOFL, index ] = TO_distance(Vstall, Clmax_to, s_ref, ... 
    Thrust, W_TO, number_engines, altitude_fi)

altitude_f = altitude_fi*0.3048;
[airDens_f, airPres_f, temp_f, soundSpeed_f] = Atmos(altitude_f);
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = ...
    convert_to_imperial(airDens_f, airPres_f, temp_f, soundSpeed_f);

[airDens_sl, airPres_sl, temp_sl, soundSpeed_sl] = Atmos(0);
[airDens_sli, airPres_sli, temp_sli, soundSpeed_sli] = ...
    convert_to_imperial(airDens_sl, airPres_sl, temp_sl, soundSpeed_sl);

sigma = airDens_fi/airDens_sli;
% Initial calculations 
Vlo = 1.2*Vstall; %kt
M_lo = Vlo * 0.00149984; %kt -> Mach
% disp(M_lo); % find thrust fraction corresponding to M_lo:
%             http://adg.stanford.edu/aa241/performance/takeoff.html

% Calculate lift of thrust. 
T_vlo = 0.75*Thrust; %Estimate on basis of thrust fit for high bypass engines.


% Main Calculations

index = (W_TO^2)/(sigma*Clmax_to*s_ref*T_vlo);

if (number_engines > 4 || number_engines < 1)
    ERROR = 'This function can not compute a take-off field length for the specified number of engines';
    disp(ERROR);
end

if (number_engines == 2)
    TOFL = 857.4 + 28.43*index + .0185*(index)^2;
end
if(number_engines == 3)
    TOFL = 667.9 + 26.91*Index + .0123*(Index)^2;
end
if(number_engines == 4)
    TOFL = 486.7 + 26.20*Index + .0093*(Index)^2;
end
end
