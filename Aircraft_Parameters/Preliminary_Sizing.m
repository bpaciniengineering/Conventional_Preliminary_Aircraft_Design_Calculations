%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Fredericks
% MAE 332 - Preliminary Sizings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Fuselage parameters from 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = .86 % From Table 6.4 Raymer, varies with type of aircraft
C = .42 % From Table 6.4 Raymer, varies with type of aircraft
Wo = 13000 % lb ---just an estimate, needs to be input from weight calcs.
fineness = 3 % length/(max diameter), 3 recommended by Raymer for subsonic aircraft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fuselage Sizing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L_fus, D_max_fus] =  Fuselage_Sizing(a, C, Wo, fineness) % ft fuselage length and max diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_W = 1; % 0/1  wings & tails geometry table  off/on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose Tail Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_Tail = 1; % 0 if a conventional tail arrangement; 1 if T-Tail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need input from outside sizings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_W =             294; % ft^2, Wing suraface area guess based on Cessna CJ3
lambda_W =        0.3; % taper ratio,per Raymer 83,assuming LE sweep is appx 30 deg
AR_W =              8; % initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimations from Raymer / Martinelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_HT =            0.9; % HT volume coefficient, Martinelli 0.9 Bizjet 1.00 transport
c_VT =           0.09; % VT volume coeffieient Initial Guess: 0.09 per Raymer pg. 160 and Martinelli
L_HT =   0.45 * L_fus; % ft, H tail moment arm : aft engine per Raymer 160: 0.45 * L_fus
L_VT =           L_HT; % ft, V tial moment arm, initial guess L_HT
AR_HT =             4; % initial guess 4 Raymer, 113
AR_VT =          0.95; % initial guess 0.95 Raymer, 113, T-Tail
lambda_HT =      0.45; % initial guess 0.45 Raymer, 113
lambda_VT =       0.8; % initial guess 0.8 Raymer, 113, T-Tail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size Wings and Tails
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W, b_HT, C_root_HT, ...
    C_tip_HT, C_bar_HT, Y_bar_HT, b_VT, C_root_VT, C_tip_VT, C_bar_VT, ...
    Y_bar_VT] = Wing_and_Tail_Sizing(S_W, lambda_W, AR_W, c_HT, ...
    L_HT, AR_HT, lambda_HT, c_VT, L_VT, AR_VT, lambda_VT, T_Tail);

Wingtype = {'Wing'; 'Horizontal Tail'; 'Vertical Tail'};
Span = [b_W; b_HT; b_VT];
Root_Chord = [C_root_W; C_root_HT; C_root_VT];
Tip_Chord = [C_tip_W; C_tip_HT; C_tip_VT];
MAC = [C_bar_W; C_bar_HT; C_bar_VT];
Y_MAC = [Y_bar_W; Y_bar_HT; Y_bar_VT];
T = table(Span, Root_Chord, Tip_Chord, MAC, Y_MAC, 'RowNames', Wingtype);
if table_W == 1
    T
end

