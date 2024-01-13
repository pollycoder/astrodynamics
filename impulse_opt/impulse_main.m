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
options = optimoptions("particleswarm", "SwarmSize", 10000, 'UseParallel', true, 'MaxIterations', 1000);
[init_X, init_result, exitflag] = particleswarm(@impulse_obj, 2, range(:, 1), range(:, 2), options);
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

% Plot
n = 1000;
timeList = linspace(0, X(2)-X(1), n)*day*365;
for i=1:length(timeList)
    [rt,~] = rv02rvf(r0,v1,timeList(i),muSun);
    r(:,i) = rt;
end

fList = linspace(0, 2*pi, n);
for i=1:length(fList)
    coeEarth = [aEarth, eEarth, iEarth, 0, 0, fList(i)];
    coeMars = [aMars, eMars, iMars, 0, 0, fList(i)];
    REarth(:,i) = coe2rv(coeEarth, muSun);
    RMars(:,i) = coe2rv(coeMars, muSun);
end

coeEarth0 = [aEarth, eEarth, iEarth, 0, 0, f0];
[rEarth0, vEarth0] = coe2rv(coeEarth0, muSun);
timeList = linspace(0, X(1), n) * day * 365;
for i=1:length(timeList)
    [rt,~] = rv02rvf(rEarth0,vEarth0,timeList(i),muSun);
    rWait(:,i) = rt;
end

plot3(REarth(1, :), REarth(2, :), REarth(3, :), 'LineWidth', 1.5, 'Color', 'b');hold on
plot3(RMars(1, :), RMars(2, :), RMars(3, :), 'LineWidth', 1.5, 'Color', 'k');hold on
plot3(r(1, :), r(2, :), r(3, :), 'r--','LineWidth', 1.5);hold on
plot3(rWait(1, :), rWait(2, :), rWait(3, :), 'm--', 'LineWidth', 1.5); hold on

plot3(rWait(1, 1), rWait(2, 1), rWait(3, 1), 'g*', 'LineWidth', 2);hold on
text(rWait(1, 1), rWait(2, 1), rWait(3, 1), 'Departure');hold on
plot3(r(1, 1), r(2, 1), r(3, 1), 'c*', 'LineWidth', 2);hold on
text(r(1, 1), r(2, 1), r(3, 1), 'Interim');hold on
plot3(r(1, end), r(2, end), r(3, end), 'r*', 'LineWidth', 2);hold on
text(r(1, end), r(2, end), r(3, end), 'Arrival');

plot3(0, 0, 0, 'k*', 'LineWidth', 3);
text(0, 0, 0, 'Sun');
axis equal

legend('Earth Orbit', 'Mars Orbit', 'Opt Trajectory', 'Waiting Orbit');
title('Trajectory Opt (PSO-Parallel)');