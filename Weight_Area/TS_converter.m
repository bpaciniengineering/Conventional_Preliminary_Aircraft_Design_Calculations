%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jan Bernhard                                                          %%
%% MAE 332 - Aircraft Design                                             %%
%% Reserved Fuel Calculations                                            %%
%% Mar. 08, 2017                                                         %%
%% Modified: xx/xx/xxxx                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE:  Converts coordinates of contraint plot into Surface area and Thrust 
%       value.

function [ S_ref, To ] = TS_converter( wing_loading, thrust_loading, W_TO )

S_ref = W_TO/wing_loading;
To = thrust_loading*W_TO;

end
