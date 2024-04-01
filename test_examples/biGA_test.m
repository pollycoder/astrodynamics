clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Big Project - Asteroid Mining
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constant
% Unit:
% Length: km
% Angle: rad
% Time: s
coeAsteroid0 = [4.374587943314110e+08, 0.134098123850821, ...
                0.0540505062211469, 2.61854482481308, ...
                4.00216803342331, 3.31673390721605];
coeEarth0 = [1.495484423703440e+08, 0.0163866660358080, ...
             5.40080930104537e-05, 3.71828887427766, ...
             4.38789863130065, 6.20499744208261];
coeMars0 = [2.279254603773820e+08, 0.0934491898618057, ...
            0.0322523881233316, 0.863747331544666, ...
            5.00261081874214, 1.94894057775148];
muMars = 4.282837521400000e+04;
muSun = 1.327124400180000e+11;
day = 86400;
g0 = 9.806650000000000e-3;
Isp = 3000;
RMars = 3.389920000000000e+03;
rpMin = 300 + RMars;

% Unit transform
% To make the calculation faster and preciser.
% From now on, all the calculation will be completed 
% in new unit system.
lUnit = 1 / coeEarth0(1);                                   % Length (AU)
tUnit = 1 / (2 * pi * sqrt(coeEarth0(1) ^ 3 / muSun));      % Time (y)
vUnit = lUnit / tUnit;                                      % Velocity (AU/y)
aUnit = lUnit / tUnit^2;                                    % Acceleration (AU/y^2)
muUnit = lUnit ^ 3 / tUnit ^ 2;                             % Mu (AU^3/y^2)
coeUnit = [lUnit, ones(1, 5)];                              % Change the unit of orbit elements quickly

% New unit - They will be set as global constant
muSunNew = muSun * muUnit;
muMarsNew = muMars * muUnit;
g0New = g0 * aUnit;
IspNew = Isp * tUnit;

coeEarth0New = coeEarth0 .* coeUnit;
coeMars0New = coeMars0 .* coeUnit;
coeAsteroid0New = coeAsteroid0 .* coeUnit;

rpMinNew = rpMin * lUnit;

tWaitUpper = 1825 * day;
tTotalUpper = 5475 * day;
tWaitUpperNew = tWaitUpper * tUnit;
lb = [1, 3, 5, 8, 11, 12, 0, 0, 0, 0]';
ub = [2, 4, 6, 9, 12, 13, 1, 1, 2 * pi, 2 * pi]';
%lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
%ub = [5, 15, 15, 15, 15, 15, 1, 1, 2 * pi, 2 * pi]';

%% Optimize fuel - global - PSO
options = optimoptions("particleswarm", "SwarmSize", 1000, ...
                       'UseParallel', true, 'MaxIterations', 1000, ...
                       'HybridFcn', 'patternsearch', 'Display', 'iter');
[X, init_result, exitflag] = particleswarm(@biGA_obj, 10, lb, ub, options);


%% Optimize fuel - local
options = optimset('MaxIter', 10000);
[X, result] = fminsearch(@biGA_obj, X, options);
fprintf("J=%f\n",result);


%% Calculate all the variables
tol = 1e-20;
penalty = 1e20;                        % Since the result should be negative, any positive number could be penalty

% Initial mass
mDry = 500;                                                 % Initial dry mass (kg)                     
mFuel = 500;                                                % Initial fuel mass (kg)
mTotal0 = mDry + mFuel;                                     % Initial total mass (kg)

% Departure: E -> M (t0-t1)
tEMNew = X(2) - X(1);                                       % Transfer time
[rE0New, vE0New] = coe2rv(coeEarth0New, muSunNew, tol);     % RV of Earth (t=0)
[rM0New, vM0New] = coe2rv(coeMars0New, muSunNew, tol);      % RV of Mars (t=0)

[rEt0New, vEt0New] = rv02rvf(rE0New, vE0New, ...
                             X(1), muSunNew);               % RV of Earth (t=t0)
[rMt1New, vMt1New] = rv02rvf(rM0New, vM0New, ...
                             X(2), muSunNew);               % RV of Mars (t=t1)

[vt0New, vt11New] = LambSol(rEt0New, rMt1New, ...
                            tEMNew, muSunNew);              % Lambert problem 1: E->M

dvt0New = vt0New - vEt0New;                                 % 1st impulse (t=t0)
dvt0NormNew = norm(dvt0New);                                % Unit transform to fit the impulse solver
dvt0Norm = dvt0NormNew / vUnit;
if dvt0Norm < 4
    dvt0New = zeros(3, 1);
    dvt0Norm = 0;
    dvt0NormNew = 0;
else
    dvt0Vector = dvt0New / dvt0NormNew;
    dvt0New = dvt0New - 4 * dvt0Vector;
    dvt0Norm = dvt0Norm - 4;
    dvt0NormNew = dvt0Norm * vUnit;
end


% Arrival: M->A (t1-t2)
tMA = X(3) - X(2);
[rA0New, vA0New] = coe2rv(coeAsteroid0New, muSunNew, tol);  % RV of Asteroid (t=0)
[rAt2New, vAt2New] = rv02rvf(rA0New, vA0New, ...
                             X(3), muSunNew);               % RV of Asteroid (t=t2)

[vt13New, vt2New] = LambSol(rMt1New, rAt2New, ...
                            tMA, muSunNew);                 % Lambert problem 2: M->A

% GA-1:SOI (t1)
[vt121New, vt122New, dvGAt1New] = SOI_opt(vt11New, vt13New, vMt1New, muMarsNew, X(7), X(9));

dvt11New = vt121New - vt11New;                              % 2nd impulse 1 (t=t1) - SOI
dvt12New = vt13New - vt122New;                              % 2nd impulse 2 (t=t1) - SOI
dvt11NormNew = norm(dvt11New);                              % Unit transform to fit the impulse solver
dvt12NormNew = norm(dvt12New);                              % Unit transform to fit the impulse solver

dvt2New = vAt2New - vt2New;                                 % 3rd impulse (t=t2)
dvt2NormNew = norm(dvt2New);                                % Unit transform to fit the impulse solver

% Return: A->M (t3-t4)
tAM = X(5) - X(4);
[rAt3New, vAt3New] = rv02rvf(rA0New, vA0New, ...
                             X(4), muSunNew);               % RV of Asteroid (t=t3)
[rMt4New, vMt4New] = rv02rvf(rM0New, vM0New, ...
                             X(5), muSunNew);               % RV of Mars (t=t4)

[vt3New, vt41New] = LambSol(rAt3New, rMt4New, ...
                            tAM, muSunNew);                 % Lambert problem 3: A->M

dvt3New = vt3New - vAt3New;                                 % 4th impulse (t=t3)
dvt3NormNew = norm(dvt3New);                                % Unit transform to fit the impulse solver

% Return: M->E (t4-t5)
tME = X(6) - X(5);
[rEt5New, vEt5New] = rv02rvf(rE0New, vE0New, ...
                             X(6), muSunNew);               % RV of Earth (t=t5)

[vt43New, vt5New] = LambSol(rMt4New, rEt5New, ...
                            tME, muSunNew);                 % Lambert problem 4: M->E

% GA-2:SOI (t4)
[vt421New, vt422New, dvGAt4New] = SOI_opt(vt41New, vt43New, vMt4New, muMarsNew, X(8), X(10));

dvt41New = vt421New - vt41New;                              % 2nd impulse 1 (t=t1) - SOI
dvt42New = vt43New - vt422New;                              % 2nd impulse 2 (t=t1) - SOI
dvt41NormNew = norm(dvt41New);                              % Unit transform to fit the impulse solver
dvt42NormNew = norm(dvt42New);                              % Unit transform to fit the impulse solver


dvt5New = vEt5New - vt5New;                                 % 6th impulse (t=t5)
dvt5NormNew = norm(dvt5New);                                % Unit transform to fit the impulse solver
dvt5Norm = dvt5NormNew / vUnit;
if dvt5Norm < 4
    dvt5New = zeros(3, 1);
    dvt5Norm = 0;
    dvt5NormNew = 0;
else
    dvt5Vector = dvt5New / dvt5NormNew;
    dvt5New = dvt5New - 4 * dvt5Vector * vUnit;
    dvt5Norm = dvt5Norm - 4;
    dvt5NormNew = dvt5Norm * vUnit;
end

% Impulses
dvt0 = dvt0New / vUnit;
dvt11 = dvt11New / vUnit;
dvt12 = dvt12New / vUnit;
dvt2 = dvt2New / vUnit;
dvt3 = dvt3New / vUnit;
dvt41 = dvt41New / vUnit;
dvt42 = dvt42New / vUnit;
dvt5 = dvt5New / vUnit;
dvt1GA = dvGAt1New / vUnit;
dvt4GA = dvGAt4New / vUnit;

% Velocity
vA0 = vA0New / vUnit;
vAt2 = vAt2New / vUnit;
vAt3 = vAt3New / vUnit;
vEt0 = vEt0New / vUnit;
vEt5 = vEt5New / vUnit;
vMt1 = vMt1New / vUnit;
vMt4 = vMt4New / vUnit;
vt0 = vt0New / vUnit;
vt11 = vt11New / vUnit;
vt121 = vt121New / vUnit;
vt122 = vt122New / vUnit;
vt13 = vt13New / vUnit;
vt2 = vt2New / vUnit;
vt3 = vt3New / vUnit;
vt41 = vt41New / vUnit;
vt421 = vt421New / vUnit;
vt422 = vt422New / vUnit;
vt43 = vt43New / vUnit;
vt5 = vt5New / vUnit;

% Rendezvous velocity - real
dvt0Real = vt0 - vEt0;
dvt5Real = vEt5 - vt5;

% Position
rA0 = rA0New / lUnit;
rAt2 = rAt2New / lUnit;
rAt3 = rAt3New / lUnit;
rE0 = rE0New / lUnit;
rEt0 = rEt0New / lUnit;
rEt5 = rEt5New / lUnit;
rM0 = rM0New / lUnit;
rMt1 = rMt1New / lUnit;
rMt4 = rMt4New / lUnit;

% Parameters
X_int = X;
X_int(1:6) = X(1:6) / tUnit / day;
X_int(7:8) = X(7:8) / lUnit;
X_int(9) = mod(X_int(9), 2 * pi);
X_int(10) = mod(X_int(10), 2 * pi);

% Mass
mSample = 1299.6120747;
[mTotalt0, dmt0] = impulseFuel(mTotal0, dvt0NormNew, IspNew, g0New);
[mTotalt11, dmt11] = impulseFuel(mTotalt0, dvt11NormNew, IspNew, g0New);
[mTotalt12, dmt12] = impulseFuel(mTotalt11, dvt12NormNew, IspNew, g0New);
[mTotalt21, dmt21] = impulseFuel(mTotalt12, dvt2NormNew, IspNew, g0New);
mTotalt22 = mTotalt21 + mSample;
[mTotalt3, dmt3] = impulseFuel(mTotalt22, dvt3NormNew, IspNew, g0New);
[mTotalt41, dmt41] = impulseFuel(mTotalt3, dvt41NormNew, IspNew, g0New);
[mTotalt42, dmt42] = impulseFuel(mTotalt41, dvt42NormNew, IspNew, g0New);
[mTotalt5, dmt5] = impulseFuel(mTotalt42, dvt5NormNew, IspNew, g0New);

mFuel_left = mFuel - dmt0- dmt11 - dmt12 - dmt21 - dmt3 - dmt41 - dmt42 - dmt5;


%% Plot
% Earth, Mars, Asteroid, Sun
plotOrbit(coeEarth0New, muSunNew);
plotOrbit(coeMars0New, muSunNew);
plotOrbit(coeAsteroid0New, muSunNew);
plot3(0,0,0,'k*','LineWidth', 3);
text(0,0,0,'Sun');

% Trajectory 1
r0 = rEt0New;
v0 = vt0New;
t01 = X(2) - X(1);
style.LineWidth = 1.5;
style.LineColor = 'r';
style.LineStyle = '--';
style.pointStyle = 'b*';
style.pointText = 'GA-Mars-1';
plotTrajectory(r0, v0, t01, muSunNew, style);
plot3(r0(1), r0(2), r0(3), 'g*', 'LineWidth', 2);
text(r0(1), r0(2), r0(3), 'Departure');

% Trajectory 2
r0 = rMt1New;
v0 = vt13New;
t12 = X(3) - X(2);
style.LineWidth = 1.5;
style.LineColor = 'm';
style.LineStyle = '--';
style.pointStyle = 'm*';
style.pointText = 'Arrival-Psyche-Mining';
plotTrajectory(r0, v0, t12, muSunNew, style);

% Trajectory 3
r0 = rAt2New;
v0 = vAt2New;
t23 = X(4)- X(3);
style.LineWidth = 1.5;
style.LineColor = 'g';
style.LineStyle = '--';
style.pointStyle = 'g*';
style.pointText = 'Return-Psyche';
plotTrajectory(r0, v0, t23, muSunNew, style);

% Trajectory 4
r0 = rAt3New;
v0 = vt3New;
t34 = X(5)- X(4);
style.LineWidth = 1.5;
style.LineColor = 'c';
style.LineStyle = '--';
style.pointStyle = 'b*';
style.pointText = 'GA-Mars-2';
plotTrajectory(r0, v0, t34, muSunNew, style);

% Trajectory 5
r0 = rMt4New;
v0 = vt43New;
t45 = X(6)- X(5);
style.LineWidth = 1.5;
style.LineColor = 'b';
style.LineStyle = '--';
style.pointStyle = 'r*';
style.pointText = 'Arrival-Earth';
plotTrajectory(r0, v0, t45, muSunNew, style);

axis equal

title("Optimized Trajectory for Psyche Mining");
legend('Earth Orbit', 'Mars Orbit', 'Psyche Orbit', '', ...
       'Earth-->Mars Trajectory', '', '', ...
       'Mars-->Psyche Trajectory', '', 'Psyche Transfer', '', ...
       'Psyche-->Mars Trajectory', '', 'Mars-->Earth Trajectory', '');
