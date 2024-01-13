%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils            
% Orbit initial valued problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt, vt] = rv02rvf(r0, v0, dt, mu)
coe = rv2coe(r0, v0, mu);
coe(6) = f0dt2ft(coe(6), dt, coe(1), coe(2), mu);
[rt, vt] = coe2rv(coe, mu);
end




