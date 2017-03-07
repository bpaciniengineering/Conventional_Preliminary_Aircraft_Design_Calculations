%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage Sizing Code - MAE 332
% Madeline Travnik
% 3 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L_fuse, D_max_fuse] = fuselage_sizing(a,C,Wo,fineness)

L_fus = a*Wo^C % ft

D_max_fuse = L_fus/fineness % ft

