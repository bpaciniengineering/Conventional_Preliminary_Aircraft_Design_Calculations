%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini & Nathan Wei                                          %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Mar. 28, 2017 Tue                                                     %%
%% Modified 07/02/2017 by Leif Fredericks (need to adjust beta product)
%% Aircraft Carpet Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  aircraft_carpetplot(M_cruise, R, AR, e, tsfc, ...
    altitude_ci, altitude_fi, loiter_dur, altitude_climbi, V_approach, ...
    V_stall, Clmax_to, Clmax_land, L_takeoff, L_landing, M_climb, ...
    rate_climb, theta_app, C_D0_c, C_DR_c, K1_c, K2_c, gamma, TR, g, ...
    carpet_x_lim, carpet_y_lim, W_TO, LD,LDC)

altitude_c = altitude_ci*0.3048;
altitude_climb = altitude_climbi*0.3048;
altitude_f = altitude_fi*0.3048;
[airDens_c, airPres_c, temp_c, soundSpeed_c] = Atmos(altitude_c); % SI
[airDens_f, airPres_f, temp_f, soundSpeed_f] = Atmos(altitude_f);
[airDens_climb, airPres_climb, temp_climb, soundSpeed_climb] = ...
    Atmos(altitude_climb);
[airDens_sl, airPres_sl, temp_sl, soundSpeed_sl] = Atmos(0);

% convert nm to miles
R = R*1.15078;

% Convert values from SI to Imperial
[airDens_ci, airPres_ci, temp_ci, soundSpeed_ci] = ...
    convert_to_imperial(airDens_c, airPres_c, temp_c, soundSpeed_c);
[airDens_fi, airPres_fi, temp_fi, soundSpeed_fi] = ...
    convert_to_imperial(airDens_f, airPres_f, temp_f, soundSpeed_f);
[airDens_climbi, airPres_climbi, temp_climbi, soundSpeed_climbi] = ...
    convert_to_imperial(airDens_climb, airPres_climb, temp_climb, ...
    soundSpeed_climb);
[airDens_sli, airPres_sli, temp_sli, soundSpeed_sli] = ...
    convert_to_imperial(airDens_sl, airPres_sl, temp_sl, soundSpeed_sl);

sigma = airDens_fi/airDens_sli;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stall
figure()
hax=axes; 
title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
hold on;
V_stall = V_stall * 1.68781; %convert to ft/s
WS_stall = ((V_stall^2)*airDens_sli*Clmax_to)/(2*g);
% Plotted later for cosmetic reasons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off
WS = linspace(1,200);
TW_takeoff = ((20.9.*WS)/(sigma*Clmax_to)).*...
    (L_takeoff-69.6.*(WS./(sigma*Clmax_to)).^(.5)).^(-1);

beta_takeoff = calculate_beta('takeoff', R, loiter_dur, 0, ...
    0, AR, e, C_D0_c+C_DR_c, tsfc, LD, LDC); 

B_TO = beta_takeoff; % Instantaneous fuel fraction thus far

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take-off - revised to include full take-off field length
index = (-28.43 + sqrt(28.43^2 - 4*(857.4-L_takeoff)*0.0185))/(2*0.0185);
% Take positive root only
% index2 = (-28.43 - sqrt(28.43^2 - 4*(857.4-L_takeoff)*0.0185))/(2*0.0185);

TW_takeoff1 = WS / (sigma*Clmax_to*index);
%TW_takeoff2 = WS / (sigma*Clmax_to*index2);
%plot(WS, TW_takeoff1, 'b');
%plot(WS, TW_takeoff2, 'c--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Cruise Flight
% TW_CCF = (beta/alpha)*(k1*(beta/q)*WS + k2 + (CD_O + CD_R)/((beta/q)*WS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb-performance
dHdt = rate_climb/60; % ft/s
V_climb = (soundSpeed_climb*3.28084)*M_climb; % ft/s
beta_climb = calculate_beta('climb', R, loiter_dur, M_cruise, ...
    soundSpeed_ci*M_cruise, AR, e, C_D0_c+C_DR_c, tsfc,LD,LDC); 
    % correct to use M_Cruise here because it is the FINAL mach #
q_climb = dynamic_pressure(airDens_climbi, V_climb, g);

alpha_climb = calculate_alpha(temp_climb, airPres_climb, temp_sl, ...
    airPres_sl, gamma, M_climb, TR);

B_climb = beta_takeoff * beta_climb; % instantaneous fuel fraction

TW_climb = (B_climb/alpha_climb).*(K1_c*(B_climb/q_climb)*WS + ...
    K2_c + (C_D0_c + C_DR_c)./((B_climb/q_climb)*WS) + ...
    (1/V_climb)*(dHdt));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cruise-performance

V_c = (soundSpeed_c*3.28084)*M_cruise; % ft/s
q_c = dynamic_pressure(airDens_ci, V_c, g);

beta_c = calculate_beta('cruise', R, loiter_dur, M_cruise, ...
    soundSpeed_ci*M_cruise, AR, e, C_D0_c+C_DR_c, tsfc,LD,LDC);

alpha_c = calculate_alpha(temp_c, airPres_c, temp_sl, airPres_sl, ...
    gamma, M_cruise, TR);

B_cruise = beta_takeoff * beta_climb * beta_c;% instantaneous fuel fraction

TW_cruise = (B_cruise/alpha_c)*(K1_c*B_cruise*WS/q_c + K2_c + ...
    (C_D0_c+C_DR_c)./(B_cruise*WS/q_c));

%n = sqrt(1 + ((V^2))/g*R); % - here we assume n = 1 (no turning)

%TW_CLT = (beta/alpha)*(k1*n^2*(beta/q)*WS + k2*n + ...
  % (CD_O + CD_R)/((beta/q)*WS));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum and minimum load factor constraints (assume at cruise)
beta_n = 1; % assume negligible fuel consumption in turn
B_turn = beta_takeoff * beta_climb * beta_c * beta_n; % inst fuel fraction
n_max = 2.1 + (24000/(W_TO+10000)); % from FAR 23.337
n_max = n_max * 1.5; % safety factor
if n_max > 3.8
    n_max = 3.8;
end
n_min = -0.4 * n_max;

TW_nmax = (B_turn/alpha_c)*(K1_c*n_max^2*WS/q_c + K2_c*n_max + ...
    (C_D0_c+C_DR_c)./(B_turn*WS/q_c));
TW_nmin = (B_turn/alpha_c)*(K1_c*n_min^2*WS/q_c + K2_c*n_min + ...
    (C_D0_c+C_DR_c)./(B_turn*WS/q_c));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loiter instantaneous fuel fraction for possible future modifications
beta_l = calculate_beta('loiter', R, loiter_dur, 0, ...
    0, AR, e, C_D0_c+C_DR_c, tsfc,LD,LDC);

B_loiter = beta_takeoff * beta_climb * beta_c * beta_n * beta_l; % unused
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landing

WS_landing = (L_landing - 50/tan(deg2rad(theta_app)))...
    *sigma*Clmax_land/79.4;
%WS_landing = ((sigma * Clmax_land)/(79.4)) *(SL  - 50/tan(theta))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot everything on carpet plot
% plot rectangles first so they're in back
rectangle('Position', [WS_landing carpet_y_lim(1) ...
    (carpet_x_lim(2)) carpet_y_lim(2)], 'FaceColor', [0 1 1]);
rectangle('Position', [WS_stall carpet_y_lim(1) ...
    carpet_x_lim(2) carpet_y_lim(2)], 'FaceColor', [1 1 0]);

% Takeoff
area(WS, TW_takeoff1, 'FaceColor', 'b');
% Stall
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 1 0]);
% Climb
area(WS, TW_climb, 'FaceColor', 'm');
% Cruise
area(WS, TW_cruise, 'FaceColor', 'g');
% Landing
line([WS_landing WS_landing], get(hax,'YLim'),'Color',[0 1 1]);
alpha(0.5); % transparency
% Load Factors
plot(WS, TW_nmax, 'r.');
plot(WS, TW_nmin, 'r--'); % thrust loading should be below this line

title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
xlim(carpet_x_lim)
ylim(carpet_y_lim)

xlim(carpet_x_lim);
ylim(carpet_y_lim);
% Takeoff
plot(WS, TW_takeoff1, 'b');
% Stall
line([WS_stall WS_stall],get(hax,'YLim'),'Color',[1 1 0]);
% Climb
plot(WS, TW_climb, 'm');
% Cruise
plot(WS, TW_cruise, 'g');
% Landing
line([WS_landing WS_landing], get(hax,'YLim'),'Color',[0 1 1]);
%alpha(0.5); % transparency
% Load Factors
plot(WS, TW_nmax, 'r-.');
plot(WS, TW_nmin, 'r--'); % thrust loading should be below this line

%plot([ws ws],  get(hax,'YLim'))


title('Constraint Plane (T/W - W/S)');
xlabel('Wing Loading [W_g/S], lb/ft^2');
ylabel('Thrust Loading [T_0/W_g]');
legend('Takeoff','Stall Velocity','Climb','Cruise','Landing','n_{max}','n_{min}');
end