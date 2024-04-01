clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An example for Bi-ellipse transfer test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muSun = 1.327e11;

aEarth = 1.496e8;
aMars = 2.279e8;
iEarth = 0;
iMars = deg2rad(1.85);

% Test
[DV, dv, dt, at, et, it] = double_ellipse(aEarth, aMars, ...
                                          iEarth, iMars, muSun);

% Plot orbits of planets
styleEarth.LineWidth = 1.5;
styleMars.LineWidth = 1.5;
styleCraft.LineWidth = 1.5;
styleCraft.PointSize = 100;
coeEarth = [aEarth, 0, iEarth, 0, 0, 0];
coeMars = [aMars, 0, iMars, 0, 0, 0];
plotOrbit(coeEarth, muSun, styleEarth); hold on
plotOrbit(coeMars, muSun, styleMars); hold on

% Plot trajectory
coeTransfer1 = [at(1), et(1), it(1), 0, 0, 0];
coeTransfer2 = [at(2), et(2), it(2), 0, 0, pi];
[r10, v10] = coe2rv(coeTransfer1, muSun);
[r20, v20] = coe2rv(coeTransfer2, muSun);
plotTrajectory(r10, v10, dt(1), muSun, false, styleCraft); hold on
plotTrajectory(r20, v20, dt(2), muSun, false, styleCraft); hold on

% Plot Sun
scatter3(0, 0, 0, 100, "red", '*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Sun');

colormap('jet')
axis equal

title('Trajectory of Bi-ellipse Transfer');
