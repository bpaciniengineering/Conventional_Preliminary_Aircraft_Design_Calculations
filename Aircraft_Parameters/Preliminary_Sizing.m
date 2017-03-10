%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Fredericks and Madeline Travnik
% MAE 332 - Preliminary Sizings 
% Dependencies: Wing_and_Tail_Sizing.m, Fuselage_Sizing.m
% Modified: 03/08/2017 by Jan Bernhard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main
clear;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rad =       0.0174533; % radian conversion
% Choose Configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_Tail =            1; % 0 if a conventional tail arrangement; 1 if T-Tail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_F =           1; % 0/1  fuselage geometry table off/on
table_W =           1; % 0/1  wings & tails geometry table  off/on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outside parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_W =             294; % ft^2, Wing suraface area guess based on Cessna CJ3
lambda_W =        0.3; % taper ratio,per Raymer 83,assuming LE sweep is appx 30 deg
AR_W =              8; % initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage parameters from Raymer and 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a =              0.67; % From Table 6.4 Raymer, varies with type of aircraft
C =              0.43; % From Table 6.4 Raymer, varies with type of aircraft
Wo =            13000; % lb ---just an estimate, needs to be input from weight calcs.
fineness =          6; % length/(max diameter), 3 recommended by Raymer for subsonic aircraft... 
                       % but 6 is more realistic for commercial 
upsweep = .436; % rad = 25 deg, max upsweep per Raymer
L_cockpit = 130/12; % ft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wing surface parameters from Raymer / Martinelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_HT =            0.9; % HT volume coefficient, Martinelli 0.9 Bizjet 1.00 transport
c_VT =           0.09; % VT volume coeffieient Initial Guess: 0.09 per Raymer pg. 160 and Martinelli
AR_HT =             4; % initial guess 4 Raymer, 113
AR_VT =          0.95; % initial guess 0.95 Raymer, 113, T-Tail
lambda_HT =      0.45; % initial guess 0.45 Raymer, 113
lambda_VT =       0.8; % initial guess 0.8 Raymer, 113, T-Tail
dihedral_W =  5 * rad; % rad, guess from Raymer Table 4.2 for low subsonic swept
dihedral_HT = 5 * rad; % deg, guess same as wing
sweep_LE_W = 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
sweep_LE_HT= 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
sweep_LE_VT= 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fuselage Sizing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L_fus, D_max_fus, L_HT, V_fus] = Fuselage_Sizing(a,C,Wo,fineness,L_cockpit,upsweep);

L_VT = L_HT; % ft, V tial moment arm, initial guess L_HT
Fuselage_Parameter = [L_fus;D_max_fus;V_fus];
fus_param_names = {'Length','Max Diameter','Volume Fuselage'};
F = table(Fuselage_Parameter, 'RowNames', fus_param_names);
if table_F == 1
    disp(F);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size Wings and Tails Initially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W, b_HT, C_root_HT, ...
    C_tip_HT, C_bar_HT, Y_bar_HT, b_VT, C_root_VT, C_tip_VT, C_bar_VT, ...
    Y_bar_VT] = Wing_and_Tail_Sizing(S_W, lambda_W, AR_W, c_HT, ...
    L_HT, AR_HT, lambda_HT, c_VT, L_VT, AR_VT, lambda_VT, T_Tail);

[sweep_c_W, sweep_TE_W] = sweep(b_W, sweep_LE_W, C_root_W, C_tip_W);
[sweep_c_HT, sweep_TE_HT] = sweep(b_HT, sweep_LE_HT, C_root_HT, C_tip_HT);
[sweep_c_VT, sweep_TE_VT] =sweep(2*b_VT, sweep_LE_VT, C_root_VT, C_tip_VT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix moment arm of Horizontail Tail and its geometries if T tail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T_Tail == 1
    dL = b_VT * tan(sweep_c_VT); % difference between L_VT and where H_VT should be
    L_HT = L_VT + dL; % set L_HT to where it should be
    [b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W, b_HT, C_root_HT, ...
    C_tip_HT, C_bar_HT, Y_bar_HT, b_VT, C_root_VT, C_tip_VT, C_bar_VT, ...
    Y_bar_VT] = Wing_and_Tail_Sizing(S_W, lambda_W, AR_W, c_HT, ...
    L_HT, AR_HT, lambda_HT, c_VT, L_VT, AR_VT, lambda_VT, T_Tail);
                      % calculate properties based on new L_HT
    [sweep_c_VT, sweep_TE_VT] =sweep(2*b_VT, sweep_LE_VT, C_root_VT, C_tip_VT);
                      % calculate sweep from new L_HT                
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find Leading Edge distance from wing root chord Leading Edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_LE_HT = Y_bar_W*tan(sweep_LE_W) + C_bar_W/4 + L_HT - [Y_bar_HT * ...
    tan(sweep_LE_HT) + C_bar_HT/4];
X_LE_VT = Y_bar_W*tan(sweep_LE_W) + C_bar_W/4 + L_VT - [Y_bar_VT * ...
    tan(sweep_LE_VT) + C_bar_VT/4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output Table for Wing Surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wingtype = {'Wing'; 'Horizontal Tail'; 'Vertical Tail'};
Span = [b_W; b_HT; b_VT];
Root_Chord = [C_root_W; C_root_HT; C_root_VT];
Tip_Chord = [C_tip_W; C_tip_HT; C_tip_VT];
MAC = [C_bar_W; C_bar_HT; C_bar_VT];
Y_MAC = [Y_bar_W; Y_bar_HT; Y_bar_VT];
LE_Sweep_deg = [sweep_LE_W / rad; sweep_LE_HT/ rad; sweep_LE_VT/ rad];
Quarter_C_Sweep_deg = [sweep_c_W/ rad; sweep_c_HT/ rad; sweep_c_VT/ rad];
TE_Sweep_deg = [sweep_TE_W/ rad; sweep_TE_HT/ rad; sweep_TE_VT/ rad];
L = [0; L_HT; L_VT]; % Length from wing 1/4 MAC to 1/4 MAC of surface
X_LE = [0; X_LE_HT; X_LE_VT]; % X position from wing root chord leading edge 
Dihedral_deg = [dihedral_W / rad; dihedral_HT / rad; 0]; 

T = table(Span, Root_Chord, Tip_Chord, MAC, Y_MAC, LE_Sweep_deg, ...
    Quarter_C_Sweep_deg,TE_Sweep_deg, L, X_LE, Dihedral_deg, 'RowNames',...
    Wingtype);
if table_W == 1
    disp(T);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Determine Trailing edge and Quarter Chord Sweeps geometrically
function [sweep_c, sweep_TE] = sweep(b, sweep_LE, C_root, C_tip)
z = b/2 * tan(sweep_LE);
x = z - (C_root / 4);
y = x + (C_tip / 4);
sweep_c = atan(y/(b/2));
sweep_TE = atan((z + C_tip - C_root) / (b/2));
end

