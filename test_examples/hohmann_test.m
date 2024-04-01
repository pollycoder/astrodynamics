clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An example for Hohmann transfer test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muSun = 1.327e11;

r1 = 1.496e8;
r2 = 2.279e8;
i1 = 0;
i2 = pi / 6;

% Test
[dv1, dv2, dv, dt, coeh] = hohmann(r1, r2, i1, i2, muSun);

% Plot orbits of planets
style1.LineWidth = 1.5;
style2.LineWidth = 1.5;
styleCraft.LineWidth = 1.5;
styleCraft.PointSize = 100;
coe1 = [r1, 0, i1, 0, 0, 0];
coe2 = [r2, 0, i2, 0, 0, 0];
plotOrbit(coe1, muSun, style1); hold on
plotOrbit(coe2, muSun, style2); hold on

% Plot trajectory
coeTransfer = [coeh(1), coeh(2), coeh(3), 0, 0, 0];
[r10, v10] = coe2rv(coeTransfer, muSun);
plotTrajectory(r10, v10, dt(1), muSun, false, styleCraft); hold on

% Plot Sun
scatter3(0, 0, 0, 100, "red", '*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Sun');

colormap('jet')
axis equal

title('Trajectory of Hohmann Transfer');
