%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% Orbit initial valued problem 
% with J2 perturbation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt, vt] = rv02rvf_aveJ2(r0, v0, dt, mu, Re)
tol = 1e-20;
J2 = 1.083e-3;                                  % J2 perturbation
coe = rv2coe(r0, v0, mu);

a = coe(1);
e = coe(2);
i = coe(3);
n = sqrt(mu / a^3);                             % Average angular velocity (s^-1)
p = a * (1 - e^2);                              % Semi-latus rectum (km)

% Average rate of change
vOmega = -3/2 * J2 * (Re / p) ^ 2 ...   
              * n * cos(i); ...
vomega = 3/4 * J2 * (Re / p) ^ 2 ...
             * n * (5 * cos(i) ^ 2 - 1);
vf = 3/4 * J2 * (Re / p) ^ 2 ...                % Caution: f = lambda + n * t, df/dt = d(lambda)/dt + n
         * sqrt(1 - e^2) * n ...
         * (2 - 3 * sin(i) ^ 2) + n;

% New orbit element
coe(4) = coe(4) + vOmega * dt;
coe(5) = coe(5) + vomega * dt;
coe(6) = coe(6) + vf * dt;

[rt, vt] = coe2rv(coe, mu, tol);
end


