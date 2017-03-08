%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage Sizing Code - MAE 332
% Madeline Travnik
% 3 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L_fus, D_max_fus, L_HT] = fuselage_sizing(a,C,Wo,fineness);

L_fus = a*Wo^C; % ft
D_max_fus = L_fus/fineness; % ft
L_HT =   0.45 * L_fus; % ft, H tail moment arm : aft engine per Raymer 160: 0.45 * L_fus

