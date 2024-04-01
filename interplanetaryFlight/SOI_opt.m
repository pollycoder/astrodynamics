%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOI: Two impulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vMid1, vMid2, dvGA] = SOI_opt(v1, v2, vPlanet, muPlanet, rp, phi)
vMid0 = [1; 1; 1];                                                                 % Use v1 as initial value
obj_func = @(vMid)SOI_impulse(vMid, v1, v2, vPlanet, muPlanet, rp, phi);
options = optimoptions("fminunc", "Display", "off");
[vMid1, ~] = fminunc(obj_func, vMid0, options);
[vMid2, dvGA] = SOI(vMid1, vPlanet, muPlanet, rp, phi);
end

% vMid: velocity right before GA
function dv = SOI_impulse(vMid, v1, v2, vPlanet, muPlanet, rp, phi)
[vMid_after, ~] = SOI(vMid, vPlanet, muPlanet, rp, phi);
dv1 = vMid - v1;
dv2 = v2 - vMid_after;
dv = norm(dv1) + norm(dv2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOI Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v2, dvGA] = SOI(v1, vPlanet, muPlanet, rp, phi)
vInf1 = v1 - vPlanet;                                                       % Residual velocity when getting in the SOI
vInf = norm(vInf1);
delta = 2 * asin(1 / (1 + vInf^2 * rp / muPlanet));

% Planet coordinate
iVector = vInf1 / norm(vInf1);
hVec = cross(v1, vPlanet);
kVector = hVec / norm(hVec);
jVector = cross(kVector, iVector);

vInf2 = vInf * (sin(delta) * cos(phi) * kVector + ...
                sin(delta) * sin(phi) * jVector + ...
                cos(delta) * iVector);

v2 = vPlanet + vInf2;
dvGA = vInf2 - vInf1;
end