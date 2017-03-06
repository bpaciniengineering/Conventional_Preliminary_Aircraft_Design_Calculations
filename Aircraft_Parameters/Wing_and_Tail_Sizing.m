%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Fredericks
% MAE 332 - Wing and Tail Geometry Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns: 
% winspan                     (b)
% root chord                  (C_root)
% tip chord                   (C_tip)
% mean aerodynamic chord      (C_bar)
% spanwise location of MAC    (Y_bar) 
% 
% for Wing (W), Horizontal Tail (HT), and Vertical Tail (VT)
%
% given:
%
% wing area                          (S_W)
% wing taper ratio                   (lambda_W)
% wing aspect ratio                  (AR_W)
% horizontal tail volume coefficient (c_HT)
% horizontal tail moment arm         (L_HT)
% horizontal tail aspect ratio       (AR_HT)
% horizontal tail taper ratio        (lambda_HT)
% horizontal tail volume coefficient (c_VT)
% horizontal tail moment arm         (L_VT)
% horizontal tail aspect ratio       (AR_VT)
% horizontal tail taper ratio        (lambda_VT)
% T-Tail configuration (0 if conventional, 1 if T-tail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W, b_HT, C_root_HT, ...
    C_tip_HT, C_bar_HT, Y_bar_HT, b_VT, C_root_VT, C_tip_VT, C_bar_VT, ...
    Y_bar_VT] = Wing_and_Tail_Sizing(S_W, lambda_W, AR_W, c_HT, ...
    L_HT, AR_HT, lambda_HT, c_VT, L_VT, AR_VT, lambda_VT, T_Tail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Areas based on wing geometries
[b_W, C_root_W, C_tip_W, C_bar_W, Y_bar_W] = ...
    wing_sizing(S_W, AR_W, lambda_W); % set geometric parameters for Wing

S_HT = SHT(c_HT, C_bar_W, S_W, L_HT, T_Tail); % set HT surface area

S_VT = SVT(c_VT, b_W, S_W, L_VT, T_Tail); % set VT surface area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric qualities of tails
[b_HT, C_root_HT, C_tip_HT, C_bar_HT, Y_bar_HT] = ...
    wing_sizing(S_HT, AR_HT, lambda_HT); % set geometric parameters for HT

[b_VT, C_root_VT, C_tip_VT, C_bar_VT,Y_bar_VT] = ...
    wing_sizing(S_VT, AR_VT, lambda_VT); % set geometric parameters for VT

Y_bar_VT = Y_bar_VT * 2; % must be doubled per Raymer 193
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% H Tail Surface area Calculation, Raymer 159
function [S_HT] = SHT(c_HT, C_bar_W, S_W, L_HT, T_Tail)
if T_Tail == 1
    c_HT = c_HT - .05 * c_HT;
end
S_HT = c_HT * C_bar_W * S_W / L_HT;
end 

% V Tail Surface area calculation, Raymer 159
function [S_VT] = SVT(c_VT, b_W, S_W, L_VT, T_Tail)
if T_Tail == 1
    c_VT = c_VT - .05 * c_VT;
end
S_VT = c_VT * b_W * S_W / L_VT;
end 

% Size wing parameters based on Raymer 192 assuming trapezoidal wing
function [b, C_root, C_tip, C_bar, Y_bar] = ...
    wing_sizing(S, AR, lambda)
b = sqrt(S * AR);
C_root = 2 * S / (b *(1+lambda));
C_tip = lambda * C_root;
C_bar = (2/3) * C_root * (1 + lambda + lambda^2) ...
    / (1 + lambda);
Y_bar = (b / 6) * (1 + 2*lambda) / (1 + lambda);
end
end
