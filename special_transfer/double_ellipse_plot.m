%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting tool for Hohmann transfer
% Input:
%   coeD: [r, i] of departure planet orbit
%   coeA: [r, i] of arrival planet orbit
%   coeH: [a, e, i] of transfer orbit
%   mu: use the sun by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function double_ellipse_plot(coeD, coeA, coet1, coet2, mu, titleName)
if nargin < 5
    titleName = '';
end
if nargin < 4
    mu = 1.327e8;
end

% Plot
f = linspace(0, 2 * pi, 1000);
ft1 = linspace(0, pi, 1000);
ft2 = linspace(-pi, 0, 1000);

% Departure and arrival
xD = zeros(3, length(f));
vD = zeros(3, length(f));
xA = zeros(3, length(f));
vA = zeros(3, length(f));
xt1 = zeros(3, length(ft1));
vt1 = zeros(3, length(ft1));
xt2 = zeros(3, length(ft2));
vt2 = zeros(3, length(ft2));
for i=1:length(f)
    coeD1 = [coeD(1), 0, coeD(2), 0, 0, f(i)];
    coeA1 = [coeA(1), 0, coeA(2), 0, 0, f(i)];
    coet11 = [coet1(1), coet1(2), coet1(3), 0, 0, ft1(i)];
    coet12 = [coet2(1), coet2(2), coet2(3), 0, 0, ft2(i)];
    [xD(:, i), vD(:, i)] = coe2rv(coeD1, mu);
    [xA(:, i), vA(:, i)] = coe2rv(coeA1, mu);
    [xt1(:, i), vt1(:, i)] = coe2rv(coet11, mu);
    [xt2(:, i), vt2(:, i)] = coe2rv(coet12, mu);
end

plot3(xD(1, :), xD(2, :), xD(3, :), 'LineWidth', 1.5, 'Color', 'b');hold on
plot3(xA(1, :), xA(2, :), xA(3, :), 'LineWidth', 1.5, 'Color', 'k');hold on
plot3(xt1(1, :), xt1(2, :), xt1(3, :), 'LineWidth', 1.2, 'Color', 'r', 'LineStyle', '-');hold on
plot3(xt2(1, :), xt2(2, :), xt2(3, :), 'LineWidth', 1.2, 'Color', 'm', 'LineStyle', '-');hold on
plot3(0, 0, 0, 'k*', 'LineWidth', 2);hold on
text(0, 0, 0, 'Sun');hold on
plot3(xt1(1, 1), xt1(2, 1), xt1(3, 1), 'g*', 'LineWidth', 2);hold on
text(xt1(1, 1), xt1(2, 1), xt1(3, 1), 'Departure');hold on
plot3(xt1(1, end), xt1(2, end), xt1(3, end), 'c*', 'LineWidth', 2);hold on
text(xt1(1, end), xt1(2, end), xt1(3, end), 'Interim');hold on
plot3(xt2(1, end), xt2(2, end), xt2(3, end), 'r*', 'LineWidth', 2);hold on
text(xt2(1, end), xt2(2, end), xt2(3, end), 'Arrival');hold on

legend('Departure Planet Orbit', 'Arrival Planet Orbit', '1st ellipse trajectory', '2nd ellipse trajectory', '', '');
title(["Double-ellipse trasfer - ", titleName]);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
end










