%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation - J2
% Input:
%   R - Position, 3x1 vector
%   mu - gravity coefficient
%   Re - radius of the planet
%%%%%%%%%%%%%%%%%%%%%%%%%
function aJ2 = perJ2(R, mu, Re)
J2 = 1.083e-3;                                  % J2 perturbation
aJ2 = zeros(3, 1);

r = norm(R);

coef1 = -3/2 * mu / r^3 * J2 * (Re / r) ^ 2;
coef2 = R(3) ^ 2 / r ^ 2;
aJ2(1) = coef1 * (1 - 5 * coef2) * R(1);
aJ2(2) = coef1 * (1 - 5 * coef2) * R(2);
aJ2(3) = coef1 * (3 - 5 * coef2) * R(3);
end