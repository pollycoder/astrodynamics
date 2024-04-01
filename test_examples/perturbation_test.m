%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation 
%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
% Constant
mu = 398600;                                    % Gravity coefficient (km^3/s^2)
Re = 6378;                                      % Radius of Earth (km)
J2 = 1.083e-3;                                  % J2 perturbation
a = 10000;                                      % Semi-major axis (km)
e = 0.1;                                        % Eccentricity
i = deg2rad(30);                                % Inclination (rad)
Omega0 = deg2rad(60);                           % Right ascension of ascending node (rad)
omega0 = deg2rad(90);                           % Argument of pericentre (rad)
f0 = deg2rad(45);                               % True anomaly (rad)

n = sqrt(mu / a^3);                             % Average angular velocity (s^-1)
p = a * (1 - e^2);                              % Semi-latus rectum (km)

% Average rate of change
vOmega = -3/2 * J2 * (Re / p) ^ 2 ...   
              * n * cos(i); ...
vomega = 3/4 * J2 * (Re / p) ^ 2 ...
             * n * (5 * cos(i) ^ 2 - 1);
vf = 3/4 * J2 * (Re / p) ^ 2 ...                % Caution: f = lambda + n * t, df/dt = d(lambda)/dt + n
         * sqrt(1 - e^2) * n ...
         * (2 - 3 * sin(i) ^ 2) + n;

% Start transform 1 - Average method
n = 1e4;
t = 240 * 3600;                                  % Transform time (24 h)

trange = linspace(0, t, n);
for j=1:length(trange)
    tj = trange(j);
    Omega = Omega0 + tj * vOmega;
    omega = omega0 + tj * vomega;
    f = f0 + tj * vf;
    coe = [a, e, i, Omega, omega, f];
    [rAve(:, j), vAve(:, j)] = coe2rv(coe, mu);
    r = rAve(:, j) ./ norm(rAve(:, j));
    v = vAve(: ,j) ./ norm(vAve(:, j));
    hAve(:, j) = cross(r, v);
end

coe(6) = mod(coe(6), 2 * pi);
rAverage = rAve(:, end);
vAverage = vAve(:, end);

coe = [a, e, i, Omega0, omega0, f0];
[r0, v0] = coe2rv(coe, mu);

% plot
figure
plotTrajectory(r0, v0, t, mu, true); hold on
title('Perturbated trajectory');


% Start transform 2 - ODE45
coe0 = [a, e, i, Omega0, omega0, f0];
[R0, V0] = coe2rv(coe0, mu);
RV0 = [R0; V0];
tspan = [0; t];
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[tArray, RV] = ode45(@(t, RV)twoBodyJ2Ode(t, RV, mu, Re), tspan, RV0, options);
rODE = RV(end, 1:3)';
vODE = RV(end, 4:6)';
coeODE = rv2coe(rODE, vODE, mu);
for j=1:size(RV, 1)
    r = RV(j, 1:3) ./ norm(RV(j, 1:3));
    v = RV(j, 4:6) ./ norm(RV(j, 4:6));
    hODE(:, j) = cross(r, v)';
end

% Output
coe
coeODE
rAverage
rODE
vAverage
vODE

%% Plot
figure
% Original orbit
f = linspace(0, 2 * pi, 1000);
for j=1:length(f)
    fj = f(j);
    coeOri = [a, e, i, Omega0, omega0, fj];
    [rOri(:, j), vOri(:, j)] = coe2rv(coeOri, mu);
    r = rOri(:, j) ./ norm(rOri(:, j));
    v = vOri(: ,j) ./ norm(vOri(:, j));
    hOri(:, j) = cross(r, v);
end

m = 1000;
% Original orbit
plot3(rOri(1, :), rOri(2, :), rOri(3, :), 'k-', 'LineWidth', 3);hold on
quiver3(zeros(size(hOri(1, 1:m:end))), zeros(size(hOri(1, 1:m:end))), ...
        zeros(size(hOri(1, 1:m:end))), hOri(1, 1:m:end), hOri(2, 1:m:end), ...
        hOri(3, 1:m:end), 10000, 'LineWidth', 2, 'Color', 'k');hold on

% Earth
plot3(0,0,0,'k*','LineWidth',5);hold on                                     

% Departure
plot3(RV0(1), RV0(2), RV0(3), 'g*', 'LineWidth', 5);hold on  

% Average method
plot3(rAve(1,end), rAve(2,end), rAve(3,end), 'r*', 'LineWidth', 5);hold on
plot3(rAve(1,:), rAve(2,:), rAve(3,:), 'r-','LineWidth', 1.5);hold on
quiver3(zeros(size(hAve(1, 1:m:end))), zeros(size(hAve(1, 1:m:end))), ...
        zeros(size(hAve(1, 1:m:end))), hAve(1, 1:m:end), hAve(2, 1:m:end), ...
        hAve(3, 1:m:end), 8000, 'LineWidth', 2, 'Color', 'm');hold on

% ODE
p = 10;
plot3(rODE(1), rODE(2), rODE(3), 'b*', 'LineWidth',5);hold on
plot3(RV(:,1),RV(:,2),RV(:,3), 'y-','LineWidth', 1.5);hold on
quiver3(zeros(size(hODE(1, 1:p:end))), zeros(size(hODE(1, 1:p:end))), ...
        zeros(size(hODE(1, 1:p:end))), hODE(1, 1:p:end), hODE(2, 1:p:end), ...
        hODE(3, 1:p:end), 5000, 'LineWidth', 2, 'Color', 'b');

legend('Original Oribit', 'h-original', 'Earth', 'Departure', ...
        'Arrival-Average', 'Average Method', ...
        'h-average', 'Arrival-ODE45', ...
        'ODE45', 'h-ODE45', 'Location', 'northeast');
title(['Comparison of the Perturbated Trajectory Generated in Average Method ' ...
        'and ode45 Solver'], 'FontSize', 15);

axis equal
grid on
