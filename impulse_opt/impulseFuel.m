%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function used to calculate fuel cost 
% during one impulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mf, dm] = impulseFuel(m0, dv, Isp, g0)
dm = m0 * (1 - exp(-(dv / (Isp * g0))));
mf = m0 - dm;
end
