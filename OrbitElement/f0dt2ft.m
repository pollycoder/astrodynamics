%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% True anomaly transform in delta t
% For parabola, replace a with p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ft = f0dt2ft(f0, dt, a, e, mu)
if nargin == 4
    mu = 398600;
end
n = sqrt(mu / a^3);                         % Mean angular velocity
if e < 0
    error("Wrong input e. e>0.");
elseif e < 1 || e > 1                       % Ellipse, circle and hyperbola
    % Initial angles
    E0 = f2E(f0, e);
    M0 = E2M(E0, e);

    % Final angles
    Mt = M0 + n * dt;
    Et = M2E(Mt, e);
    ft = E2f(Et, e);
else                                        % Parabola
    % Barker's equation
    B = @(t)3 * n * t;                      % Find perigee moment
    obj_func = @(t)(B(t) + sqrt(1 + B(t)^2)) ^ (1 / 3)+ ...
                  (B(t) - sqrt(1 + B(t)^2)) ^ (1 / 3) - tan(0.5 * f0);
    t0 = newton(obj_func, 0);

    t = t0 + dt;
    barker = @(t)(B(t) + sqrt(1 + B(t)^2)) ^ (1 / 3)+ ...
               (B(t) - sqrt(1 + B(t)^2)) ^ (1 / 3);
    ft = 2 * atan(barker(t));
end

end