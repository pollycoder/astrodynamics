%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Double-ellipse bi-impulse transfer
% Treat the orbit as circular orbit.
% Inclination can be nonzero.
% CAUTION:
%   Only suitable for interplanetary flight !!!
% Input:
%   r1: radius of the departure planet 
%   r2: radius for the arrival planet 
%   i1: inclination of the departure planet
%   i2: inclination of the arrival planet
%   mu: gravity coefficient of the central body, 
%       use the sun by default (1.327e11 km3/s2)
% Output:
%   DV(1): departure impulse 
%   DV(2): interim impulse
%   DV(3): arrival impulse
%   dv: total impulse
%   dt: transfer time
%   at(1): semi-major axis of the 1st transfer orbit
%   et(1): eccentricity of the 1st trasfer orbit
%   it(1): inclination of the 1st trasfer orbit
%   at(2): semi-major axis of the 2nd transfer orbit
%   et(2): eccentricity of the 2nd trasfer orbit
%   it(2): inclination of the 2nd trasfer orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DV, dv, dt, at, et, it] = double_ellipse(r1, r2, i1, i2, mu)
if nargin < 5
    mu = 1.327e11;
end

% Velocities of the two planets
v1 = sqrt(mu / r1);
v2 = sqrt(mu / r2);

% Optimizing variable: at1, it1, it2
syms a1 ih1 ih2
at2 = a1 + abs(r2 - r1) / 2;
ri = 2 * a1 - r1;
Et1 = -mu / (2 * a1);
Et2 = -mu / (2 * at2);

vt11 = sqrt(2 * (Et1 + mu / r1));
vt12 = sqrt(2 * (Et1 + mu / ri));
vt21 = sqrt(2 * (Et2 + mu / ri));
vt22 = sqrt(2 * (Et2 + mu / r2));

% Impulse
dv1 = sqrt(v1^2 + vt11^2 - 2 * v1 * vt11 * cos(ih1 - i1));
dv2 = sqrt(vt12^2 + vt21^2 - 2 * vt12 * vt21 * cos(ih1 - ih2));
dv3 = sqrt(vt22^2 + v2^2 - 2 * vt22 * v2 * cos(ih2 - i2));
dv = dv1 + dv2 + dv3;

% 1st order necessary condition
%d2v = [diff(dv, a1); diff(dv, ih1); diff(dv, ih2)];
dv = @(x)double(subs(dv, {a1, ih1, ih2}, {x(1), x(2), x(3)}));
x0 = [2 * (r1 + r2); i1; i2];
output = fminunc(dv, x0);

% Transfer orbit 1
at(1) = output(1);
it(1) = output(2);
et(1) = (at(1) - r1) / at(1);
T1 = 2 * pi / sqrt(mu / at(1)^3);

% Transfer orbit 2
at(2) = double(subs(at2, a1, at(1)));
it(2) = output(3);
et(2) = (at(2) - r2) / at(2);
T2 = 2 * pi / sqrt(mu / at(2)^3);

% Impulse
DV(1) = double(subs(dv1, {a1, ih1, ih2}, {at(1), it(1), it(2)}));
DV(2) = double(subs(dv2, {a1, ih1, ih2}, {at(1), it(1), it(2)}));
DV(3) = double(subs(dv3, {a1, ih1, ih2}, {at(1), it(1), it(2)}));
dv = sum(DV);
if ~isreal(dv)
    error("Double-ellipse transfer is not suitable for this flight.");
end

% Transfer time
dt = (T1 + T2) / 2;

end