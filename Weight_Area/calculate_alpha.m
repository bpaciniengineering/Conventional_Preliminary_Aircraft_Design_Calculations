%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nathan Wei                                                            %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% 7 March 2017                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = calculate_alpha(temp, airPres, temp_sl, ...
    airPres_sl, gamma, M, TR)

% calculate alpha_tilde (for high bypass ratio turbofan engine) - MOVE
theta0 = (temp/temp_sl)*(1+0.5*(gamma-1)*M^2);
delta0 = (airPres/airPres_sl)*...
    (1+0.5*(gamma-1)*M^2)^(gamma/(gamma-1));

if theta0 <= TR
    alpha = delta0*(1 - 0.49*M^0.5);
else
    alpha = delta0*(1 - 0.49*M^0.5 - 3*(theta0-TR)/(1.5+M));
end

end

