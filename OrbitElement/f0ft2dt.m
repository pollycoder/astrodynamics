%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% True anomaly transform in delta t
% For parabola, replace a with p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = f0ft2dt(f0, ft, a, e, mu)
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
    Et = f2E(ft, e);
    Mt = E2M(Et, e);

    dt = (Mt - M0) / n;
else                                        % Parabola
    % Barker's equation
    B = @(f)tan(0.5 * f) * (tan(0.5 * f) ^ 2 + 3);
    B0 = B(f0);
    Bt = B(ft);
    dt = sqrt(2) / 3 * (Bt - B0) / n;
end

end