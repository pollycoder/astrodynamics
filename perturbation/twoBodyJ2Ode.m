%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-body motion equation with J2 Perturbation
% Constant: mu, Re, J2.
% To make it easier to unify the unit, 
% set the 3 constants as inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRVdt = twoBodyJ2Ode(t, RV, mu, Re)
dRVdt = zeros(6, 1);
dRVdt(1:3) = RV(4:6);                           % dr = v

% Perturbation acceleration
aJ2 = perJ2(RV(1:3), mu, Re);

% Two-body equation with perturbation
dRVdt(4:6) = -mu / r^3 * RV(1:3) + aJ2;
end