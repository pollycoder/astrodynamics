%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% Mean anomaly ---> Eccentric anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = M2E(M, e)
if e < 0
    error("Wrong input e. e>=0.");
elseif e < 1                                    % Ellipse and circle
    kepler = @(E)E - M - e * sin(E);
    E = newton(kepler, M);                      % Newton iteration method
elseif e == 1                                   % Parabola
    fprintf("Parabola. E=M.");
    E = M;
else                                            % Hyperbola
    kepler = @(E)E + M - e * sinh(E);
    E = newton(kepler, asinh(M / e));
end

end