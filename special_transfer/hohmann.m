%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hohmann Transfer
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
%   dv1: departure impulse 
%   dv2: arrival impulse
%   dv: total impulse
%   dt: transfer time
%   ah: semi-major axis of the transfer orbit
%   eh: eccentricity of the trasfer orbit
%   ih: inclination of the trasfer orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dv1, dv2, dv, dt, ah, eh, ih] = hohmann(r1, r2, i1, i2, mu)
if nargin < 5
    mu = 1.327e11;
end

% Velocities of the two planets
v1 = sqrt(mu / r1);
v2 = sqrt(mu / r2);

% Trasfer orbit
ah = (r1 + r2) / 2;
eh = abs(r1 - r2) / (r1 + r2);
Eh = -mu / (2 * ah);
v1h = sqrt(2 * (Eh + mu / r1));
v2h = sqrt(2 * (Eh + mu / r2));

% Transfer time
n = sqrt(mu / ah^3);
T = 2 * pi / n;
dt = T / 2;

% Inclination of the transfer orbit
% Use 1st order necessary condition to find the optimal inclination
syms theta

% Impulse
dv1 = sqrt(v1^2 + v1h^2 - 2 * v1 * v1h * cos(theta - i1));
dv2 = sqrt(v2^2 + v2h^2 - 2 * v2 * v2h * cos(theta - i2));
dv = dv1 + dv2;

% 1st order condition => inclination
d2v = diff(dv, theta);
d2v = @(i)subs(d2v, theta, i);
ih = newton(d2v, i1, 1e-12, 1000);

% Optimal impulse
dv1 = double(subs(dv1, theta, ih));
dv2 = double(subs(dv2, theta, ih));
dv = double(subs(dv, theta, ih));

end