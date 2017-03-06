%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jan Bernhard                                                          %%
%% MAE 332 - Aircraft Design                                             %%
%% Reserved Fuel Calculations                                            %%
%% Feb. 23, 2017 Thur                                                    %%
%% Modified: xx/xx/xxxx                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_reservedFuel] = fuelReserve(WTO,AR,tsfc,rangeCondition,M_cruise,...
    altCruise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cruise speed
[airDens, airPres, temp, soundSpeed] = Atmos(altCruise); %kg/m^3 N/m^2 K m/s
soundSpeed = soundSpeed*2.23694;    %meter/sec -> mph
V_cruise = M_cruise*(soundSpeed);   %miles/hr
R = rangeCondition*1.15078;         %nm -> m
loiter_dur = .75;                   %hours


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Failed landing approach.
ratio_landing_1 = 0.9725;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start through
ratio_startup = 0.9725;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise climb
if M_cruise < 1
    ratio_climb = 1.0065 - 0.0325*M_cruise;                     %ESTIMATE
else
    ratio_climb = 0.991 - 0.007*M_cruise - 0.01*M_cruise^2;
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
ratio_CO = 1/(exp(R*((tsfc)/(V_cruise))/(L_D)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loiter
E = loiter_dur/(3600);
L_D_max = L_D / 0.94;
ratio_Loiter = 1/(exp((E * tsfc)/(L_D_max)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing
ratio_landing_2 = 0.9725; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Ratio
weight_ratio = ratio_landing_1 * ratio_startup * ratio_climb * ratio_CO...
    * ratio_Loiter * ratio_landing_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_reservedFuel = (WTO - (WTO*weight_ratio));
end
