%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-body motion equation
% Constant: mu, Re, J2.
% To make it easier to unify the unit, 
% set the 3 constants as inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRVdt = twoBodyOde(t, RV, mu)
dRVdt = zeros(6, 1);
dRVdt(1:3) = RV(4:6);                           % dr = v
r = norm(RV(1:3));
dRVdt(4:6) = -mu / r^3 * RV(1:3);
end