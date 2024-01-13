%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting tool for Hohmann transfer
% Input:
%   coeD: [r, i] of departure planet orbit
%   coeA: [r, i] of arrival planet orbit
%   coeH: [a, e, i] of transfer orbit
%   mu: use the sun by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hohmann_plot(coeD, coeA, coeH, mu, titleName)
if nargin < 5
    titleName = '';
end
if nargin < 4
    mu = 1.327e8;
end
% Plot
f = linspace(0, 2 * pi, 10000);
ft = linspace(0, pi, 10000);

% Departure and arrival
xD = zeros(3, length(f));
vD = zeros(3, length(f));
xA = zeros(3, length(f));
vA = zeros(3, length(f));
xH = zeros(3, length(f));
vH = zeros(3, length(f));
for i=1:length(f)
    coeD1 = [coeD(1), 0, coeD(2), 0, 0, f(i)];
    coeA1 = [coeA(1), 0, coeA(2), 0, 0, f(i)];
    coeH1 = [coeH(1), coeH(2), coeH(3), 0, 0, ft(i)];
    [xD(:, i), vD(:, i)] = coe2rv(coeD1, mu);
    [xA(:, i), vA(:, i)] = coe2rv(coeA1, mu);
    [xH(:, i), vH(:, i)] = coe2rv(coeH1, mu);
end

plot3(xD(1, :), xD(2, :), xD(3, :), 'LineWidth', 1.5, 'Color', 'b');hold on
plot3(xA(1, :), xA(2, :), xA(3, :), 'LineWidth', 1.5, 'Color', 'k');hold on
plot3(xH(1, :), xH(2, :), xH(3, :), 'LineWidth', 1.2, 'Color', 'r', 'LineStyle', '-');hold on
plot3(0, 0, 0, 'k*', 'LineWidth', 2);hold on
text(0, 0, 0, 'Sun');hold on
plot3(xH(1, 1), xH(2, 1), xH(3, 1), 'g*', 'LineWidth', 2);hold on
text(xH(1, 1), xH(2, 1), xH(3, 1), 'Departure');hold on
plot3(xH(1, end), xH(2, end), xH(3, end), 'r*', 'LineWidth', 2);hold on
text(xH(1, end), xH(2, end), xH(3, end), 'Arrival');
legend('Departure Planet Orbit', 'Arrival Planet Orbit', 'Hohmann transfer trajectory', '', '', '');
title(["Hohmann trasfer - ", titleName]);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
end










