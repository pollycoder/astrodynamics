%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOI Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v2, dvGA] = SOI_after(v1, vPlanet, muPlanet, rp, phi)
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