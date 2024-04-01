%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the trajectory from r0, v0
% The length unit should be kilometer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTrajectory(r0, v0, dt, mu, ifper, style)
if nargin == 5
    sz = 8;
    width = 1.5;
elseif nargin == 4
    sz = 8;
    width = 1.5;
    ifper = false;
elseif nargin == 6
    sz = style.PointSize;
    width = style.LineWidth;
else
    error('Not enough inputs.');
end

if ~ifper
    n = 1e3;
    t = linspace(0, dt, n);
else
    n = 1e5;
    t = linspace(0, dt, n);
end

for i=1:length(t)
    if ~ifper
        [r(:, i), ~] = rv02rvf(r0, v0, t(i), mu);
    else
        Re = 6378;
        [r(:, i), ~] = rv02rvf_aveJ2(r0, v0, t(i), mu, Re);
    end
end

rNorm = vecnorm(r);
color = r(1, :) ./ rNorm;

% Trajectory
scatter3(r(1,:), r(2,:), r(3,:), width, color, 'filled');hold on

% Departure and arrival
scatter3(r(1, 1), r(2, 1), r(3, 1), sz, color(1), 'filled');hold on
scatter3(r(1, end), r(2, end), r(3, end), sz, color(end), 'filled');hold on
end