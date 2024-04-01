%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% Eccentric anomaly ---> Mean anomaly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = E2M(E, e)
if e < 0
    error("Wrong input e. e>=0.");
elseif e < 1
    M = E - e * sin(E);
elseif e == 1
    fprintf("Parabola. Let M=E.");
    M = E;
else
    M = e * sinh(E) - E;
end

end