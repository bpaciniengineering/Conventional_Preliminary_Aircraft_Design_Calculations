%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernardo Pacini                                                       %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% Feb. 20, 2017 Mon                                                     %%
%%                                                                       %%
%% Description: This code takes in an altitude in meters and outputs the %%
%% dynamic viscosity of air at that altitude                             %%
%%                                                                       %%
%% Source:                                                               %%
%%    http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q] = dynamic_viscosity(alt)
%altitude [m]
Z = [-1000 0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 ...
    20000 25000 30000 40000 50000 60000 70000 80000];
% Dynamic Viscosity (10^-5) [N s/m^2]
mu = [1.821 1.789 1.758 1.726 1.694 1.661 1.628 1.595 1.561 1.527 1.493 ...
    1.458 1.422 1.422 1.448 1.475 1.601 1.704 1.584 1.438 1.321]*1e-5;
q = interp1(Z, mu, alt);
end
