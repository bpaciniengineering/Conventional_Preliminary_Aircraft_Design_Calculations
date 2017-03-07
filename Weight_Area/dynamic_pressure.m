%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nathan Wei                                                            %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% 7 March 2017                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = dynamic_pressure(airDens, V, g)

% Imperial units dynamic pressure calculation (density in lbm/ft^3)
if g > 30
    q = 0.5*(airDens/g)*V^2; % note: 1 lbm = 1/g slugs
% Metric units calculation (density in kg/m^3, no scaling by g necessary)
else
    q = 0.5*airDens*V^2;
end

end

