%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage Sizing Code - MAE 332
% Madeline Travnik
% 3 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L_fus, D_max_fus] = fuselage_sizing(a,C,Wo,fineness);

L_fus = a*Wo^C; % ft

D_max_fus = L_fus/fineness; % ft

