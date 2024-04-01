%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% True anomaly ---> Eccentric anomaly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = f2E(f, e)
if e < 0
    error("Wrong input e. e>=0.");
elseif e < 1                                % Ellipse and circle
    f0 = mod(f, 2*pi);
    if f0 > pi
        f0 = f0 - 2 * pi;
    end
    if f0 < -pi
        f0 = f0 + 2 * pi;
    end
    E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(0.5 * f));
    E = E + f - f0;
elseif e == 1                               % Parabola
    fprintf("Parabola. E=f");
    E = f;
else                                        % Hyperbola
    % Pay attention that the f range of 
    % hyperbola is smaller than [0,2*pi)
    ftemp = mod(f, 2 * pi);
    bound = acos(1 / e);
    if (ftemp < pi && ftemp > bound) || ...
       (ftemp > pi && 2 * pi - ftemp > bound)
        error("Wrong input f. The hyperbola orbit is impossible.");
    else
        E = 2 * atanh(sqrt((e - 1) / (e + 1)) * tanh(0.5 * f));
    end   
end

end