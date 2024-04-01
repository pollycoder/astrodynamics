%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the orbit from orbit elements - whole loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotOrbit(coe0, mu, style)
if nargin == 2
    width = 1.5;
    %color = [1, 10];
elseif nargin == 3
    width = style.LineWidth;
    %color = style.ColorScale;
else
    error('Not enough inputs.');
end

n = 1e4;
tol = 1e-20;
f = linspace(0, 2 * pi, n);
for i = 1:length(f)
    coe = coe0;
    coe(6) = f(i);
    r(:, i) = coe2rv(coe, mu, tol);
end

plotTrajectory_r(r, width);
end