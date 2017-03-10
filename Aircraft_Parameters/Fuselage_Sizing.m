%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuselage Sizing Code - MAE 332                 %
% Madeline Travnik and Jesús Serrano Cendejas    %
% 3 March 2017                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Requires a, C, Wo, fineness ratio, length of cockpit, and upsweep angle.
% Returns fuselage length, maximum fuselage diameter, moment arm of horizontal tail,
% and volume of fuselage (using top and side area estimations).

function [L_fus, D_max_fus, L_HT, V_fus] = ...
                Fuselage_Sizing(a,C,Wo,fineness,L_cockpit,upsweep)

L_fus = a*Wo^C; % ft
D_max_fus = L_fus/fineness; % ft
L_HT =   0.45 * L_fus; % ft, H tail moment arm : aft engine per Raymer 160

% Top View Calculation
A_top = D_max_fus*L_fus-2*(.5 * (.4 * D_max_fus)*L_cockpit); 
                     % ft, assuming that the diameter of the nose is 1/5 
                     % diameter of fuselage from top view

% Side View Calculation
L_fus_angled = D_max_fus/ tan(upsweep); % ft, horizontal distance from rear
                                       % where fuselage angles upward
A_side = D_max_fus*L_fus-(.5*D_max_fus*L_fus_angled)-...
    2*(.5*(.6*D_max_fus)*L_cockpit);
                     % ft, assuming that the diameter of the nose is 3/5 
                     % diameter of fuselage from side view

                     
% Volume Calculation
V_fus = 3.4 * A_top * A_side / (4*L_fus);
end

