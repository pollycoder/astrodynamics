%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative motion equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRVdt = relativeODE(t, RV)
omega = 1.135e-3;                       % Angular velocity of chief
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
A = [zeros(3), eye(3); M1, M2];
dRVdt = A * RV;
end