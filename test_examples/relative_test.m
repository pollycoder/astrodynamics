clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test example for relative motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = [0, 20000];
n = 1.135e-3;
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);

% Periodic solution: C4 = 2n*x(t0) + (dy/dt)|t0 = 0
X0_p = -165.5;
Y0_p = 672.8;
Z0_p = 10;
V_x0_p = 0.646;
V_y0_p = -2 * n * X0_p;
V_z0_p = 0.89;
RV0_period = [X0_p; Y0_p; Z0_p; V_x0_p; V_y0_p; V_z0_p];
[t_period, RV_period] = ode45(@relativeODE, tspan, RV0_period, options);
RV_period = RV_period';

C1_p = -(3 * X0_p + 2 * V_y0_p / n);
C2_p = V_x0_p / n;
C3_p = Y0_p - 2 * V_x0_p / n;
C4_p = 2 * n * X0_p + V_y0_p;
C5_p = Z0_p;
C6_p = V_z0_p / n;

fprintf('*******************************************')
fprintf('*********** Periodic Solution *************')
fprintf('*******************************************')
C1_p
C2_p
C3_p
C4_p
C5_p
C6_p


% Aperiodic solution
pert = 1e-1;
X0_a = -165.5;
Y0_a = 672.8;
Z0_a = 10;
V_x0_a = 0.646;
V_y0_a = -2 * n * X0_a + pert;
V_z0_a = 0.89;
RV0_aperiod = [X0_a; Y0_a; Z0_a; V_x0_a; V_y0_a; V_z0_a];
[t_aperiod, RV_aperiod] = ode45(@relativeODE, tspan, RV0_aperiod, options);
RV_aperiod = RV_aperiod';

C1_a = -(3 * X0_a + 2 * V_y0_a / n);
C2_a = V_x0_a / n;
C3_a = Y0_a - 2 * V_x0_a / n;
C4_a = 2 * n * X0_a + V_y0_a;
C5_a = Z0_a;
C6_a = V_z0_a / n;

fprintf('*******************************************')
fprintf('********** Aperiodic Solution *************')
fprintf('*******************************************')
C1_a
C2_a
C3_a
C4_a
C5_a
C6_a

% Plot
% Final trajectory
figure
plot3(RV_period(1, :), RV_period(2, :), RV_period(3, :), 'LineWidth', 1.5); hold on
plot3(RV_aperiod(1, :), RV_aperiod(2, :), RV_aperiod(3, :), 'LineWidth', 1.5);
title('Trajectories of periodic and aperiodic solution')
legend('Periodic', 'Aperiodic');
axis equal
grid on

% X Y Z plot
figure
subplot(1, 3, 1)
plot(t_period, RV_period(1, :), 'LineWidth', 1.5);hold on
plot(t_aperiod, RV_aperiod(1, :), 'LineWidth', 1.5);
legend('Periodic', 'Aperiodic');
xlabel('t')
ylabel('x')
title('X-t plot');

subplot(1, 3, 2)
plot(t_period, RV_period(2, :), 'LineWidth', 1.5);hold on
plot(t_aperiod, RV_aperiod(2, :), 'LineWidth', 1.5);
legend('Periodic', 'Aperiodic');
xlabel('t')
ylabel('y')
title('Y-t plot');

subplot(1, 3, 3)
plot(t_period, RV_period(3, :), 'LineWidth', 1.5);hold on
plot(t_aperiod, RV_aperiod(3, :), 'LineWidth', 1.5);
legend('Periodic', 'Aperiodic');
xlabel('t')
ylabel('z')
title('Z-t plot');

sgtitle('X Y Z transformation of periodic and aperiodic trajectories')
