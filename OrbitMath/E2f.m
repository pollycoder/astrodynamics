%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% Eccentric anomaly ---> True anomaly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = E2f(E, e)
if e < 0                                        % Wrong input
    error("Wrong eccentricity. e>0.");
elseif e < 1                                    % Circle and ellipse
    Etemp = mod(E, 2 * pi);                     % E range: [-pi, +pi]
    if Etemp > pi
        Etemp = Etemp - 2 * pi;
    end
    if Etemp < -pi
        Etemp = Etemp + 2 * pi;
    end
    f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(0.5 * Etemp));
    f = f + E - Etemp;
elseif e > 1                                    % Hyperbola
    f = 2 * atan(sqrt((1 + e) / (e - 1)) * tanh(0.5 * E));
else
    f = E;
    fprintf("Parabola, no eccentric anomaly. Let f=E.");
end

end
