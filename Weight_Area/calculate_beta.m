%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nathan Wei                                                            %%
%% MAE 332 - Aircraft Design                                             %%
%% Preliminary Design Calculations                                       %%
%% 7 March 2017                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R (miles), E (hrs), V (mph), tsfc (lbfuel/hr/lbthrust)
function beta = calculate_beta(ID, R, E, M, V, AR, e, C_D0, tsfc)


if M < 1
    L_D = AR + 10;                                              %M<1
else
    L_D = 11/sqrt(M);                                           %M>1
end

switch ID
    
    case 'takeoff'
        beta = 0.9725; % estimate
        
    case 'climb'
        if M < 1
            beta = 1.0065 - 0.0325*M;
        else
            beta = 0.991 - 0.007*M - 0.01*M^2;
        end
        
    case 'cruise'
        % Breguet Range Equation
        % R = (V/tsfc) * (L_D) * ln(Wi/Wf) %lbfuel/h/lbt
        %disp(R)
        %disp(tsfc)
        %disp(V)
        %disp(L_D)
        beta = 1/(exp(R*((tsfc/V))/(L_D)));
        
    case 'loiter'
        %LD_max = 0.5*sqrt(pi*e*AR/C_D0);
        LD_max = L_D / 0.94;
        beta = 1/(exp((E*tsfc)/LD_max));
        
    case 'land'
        beta = 0.9725;
        
end

end

