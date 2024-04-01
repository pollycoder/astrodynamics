clc
clear all
%%%%%%%%%%%%%
%% Constant
%%%%%%%%%%%%%
day = 86400;
muSun = 1.327e11;
aEarth = 1.496e8;
aMars = 2.279e8;
iEarth = 0;
iMars = deg2rad(1.85);
eEarth = 0.0167;
eMars = 0.0549;

%%%%%%%%%%%%%%%%%%%
% Two-body problem
%%%%%%%%%%%%%%%%%%%
range=[0,5;0,5];
options = optimoptions("particleswarm", "SwarmSize", 1000, ...
                       'UseParallel', true, 'MaxIterations', 1000);
[init_X, init_result, exitflag] = particleswarm(@impulse_obj, 2, ...
                                                range(:, 1), range(:, 2), ...
                                                options);
[X, result] = fminsearch(@impulse_obj, init_X);

%%%%%%%%%%%%%%%
%% Plot
%%%%%%%%%%%%%%%
% Set the time
delta_X=X(2)-X(1);                      % Normalized transfer_time
wait_time=X(1)*day*365;                 % Waiting time (s)
transfer_time=delta_X*day*365;          % Transfer time (s)
total_time=X(2)*day*365;                % Total time (s)

f0 = 0;
ft = 2/ 3 * pi;
f0t = f0dt2ft(f0, wait_time, aEarth, 0, muSun);
coeEarth = [aEarth, eEarth, iEarth, 0, 0, f0t];
ftt = f0dt2ft(ft, total_time, aMars, 0, muSun);
coeMars = [aMars, eMars, iMars,0, 0, ftt];

% Initial position
[r0,v0]=coe2rv(coeEarth, muSun);

% Final position
[rf,vf]=coe2rv(coeMars, muSun);

% Lambert
[v1,v2,~,~,~,~]=LambSol(r0,rf,transfer_time,muSun);


%% Plot
clc

% Style of each trajectory
styleEarth.LineWidth = 1.5;
styleMars.LineWidth = 1.5;
styleCraft.LineWidth = 1.5;
styleCraft.PointSize = 100;

plotOrbit(coeEarth, muSun, styleEarth);hold on
plotOrbit(coeMars, muSun, styleMars);hold on
plotTrajectory(r0, v1, transfer_time, muSun, styleCraft);hold on
colormap('jet');

scatter3(0, 0, 0, 100, "red", '*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Sun');
axis equal

title('Trajectory Opt (PSO-Parallel)');