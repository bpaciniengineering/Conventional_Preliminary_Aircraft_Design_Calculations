%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Fredericks and Madeline Travnik
% MAE 332 - Preliminary Sizings 
% Dependencies: Wing_and_Tail_Sizing.m, Fuselage_Sizing.m, 
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
% Fuselage parameters from 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a =              0.67; % From Table 6.4 Raymer, varies with type of aircraft
C =              0.43; % From Table 6.4 Raymer, varies with type of aircraft
Wo =            13000; % lb ---just an estimate, needs to be input from weight calcs.
fineness =          3; % length/(max diameter), 3 recommended by Raymer for subsonic aircraft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wing surface parameters from Raymer / Martinelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_HT =            0.9; % HT volume coefficient, Martinelli 0.9 Bizjet 1.00 transport
c_VT =           0.09; % VT volume coeffieient Initial Guess: 0.09 per Raymer pg. 160 and Martinelli
L_HT =   0.45 * L_fus; % ft, H tail moment arm : aft engine per Raymer 160: 0.45 * L_fus
L_VT =           L_HT; % ft, V tial moment arm, initial guess L_HT
AR_HT =             4; % initial guess 4 Raymer, 113
AR_VT =          0.95; % initial guess 0.95 Raymer, 113, T-Tail
lambda_HT =      0.45; % initial guess 0.45 Raymer, 113
lambda_VT =       0.8; % initial guess 0.8 Raymer, 113, T-Tail
diherdral_W = 5 * rad; % rad, guess from Raymer Table 4.2 for low subsonic swept
dihedral_HT = 0 * rad; % deg, keep zero for now
sweep_LE_W = 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
sweep_LE_HT= 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
sweep_LE_VT= 30 * rad; % rad, guess assuming M=8 from Raymer fig 4.19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fuselage Sizing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L_fus, D_max_fus] =  Fuselage_Sizing(a, C, Wo, fineness); % ft fuselage length and max diameter

Fuselage_Parameter = [L_fus;D_max_fus];
fus_param_names = {'Length','Max Diameter'};
F = table(Fuselage_Parameter, 'RowNames', fus_param_names);
if table_F == 1
    F
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size Wings and Tails
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W, b_HT, C_root_HT, ...
    C_tip_HT, C_bar_HT, Y_bar_HT, b_VT, C_root_VT, C_tip_VT, C_bar_VT, ...
    Y_bar_VT] = Wing_and_Tail_Sizing(S_W, lambda_W, AR_W, c_HT, ...
    L_HT, AR_HT, lambda_HT, c_VT, L_VT, AR_VT, lambda_VT, T_Tail);

[sweep_c_W, sweep_TE_W] = sweep(b_W, sweep_LE_W, C_root_W, C_tip_W);
[sweep_c_HT, sweep_TE_HT] = sweep(b_HT, sweep_LE_HT, C_root_HT, C_tip_HT);
[sweep_c_VT, sweep_TE_VT] =sweep(2*b_VT, sweep_LE_VT, C_root_VT, C_tip_VT);

Wingtype = {'Wing'; 'Horizontal Tail'; 'Vertical Tail'};
Span = [b_W; b_HT; b_VT];
Root_Chord = [C_root_W; C_root_HT; C_root_VT];
Tip_Chord = [C_tip_W; C_tip_HT; C_tip_VT];
MAC = [C_bar_W; C_bar_HT; C_bar_VT];
Y_MAC = [Y_bar_W; Y_bar_HT; Y_bar_VT];
LE_Sweep_deg = [sweep_LE_W / rad; sweep_LE_HT/ rad; sweep_LE_VT/ rad];
Quarter_C_Sweep_deg = [sweep_c_W/ rad; sweep_c_HT/ rad; sweep_c_VT/ rad];
TE_Sweep_deg = [sweep_TE_W/ rad; sweep_TE_HT/ rad; sweep_TE_VT/ rad];
Dihedral

T = table(Span, Root_Chord, Tip_Chord, MAC, Y_MAC, LE_Sweep_deg, ...
    Quarter_C_Sweep_deg,TE_Sweep_deg, 'RowNames', Wingtype);
if table_W == 1
    T
end




function [sweep_c, sweep_TE] = sweep(b, sweep_LE, C_root, C_tip)
z = b/2 * tan(sweep_LE);
x = z - (C_root / 4);
y = x + (C_tip / 4);
sweep_c = atan(y/(b/2));
sweep_TE = atan((z + C_tip - C_root) / (b/2));
end
