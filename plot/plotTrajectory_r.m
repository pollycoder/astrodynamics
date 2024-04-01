%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot trajectory based on ODE solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTrajectory_r(r, width)
rNorm = vecnorm(r);
color = r(1, :) ./ rNorm;
scatter3(r(1,:), r(2,:), r(3,:), width, color, 'filled');
end