% MEGN505
% Final Project
% John Blessinger & Adrian Davis
% crash_test_dynamics_model.m

%% MAIN

clc; clear all; close all;

% Input parameters

% Body dimensions
A = 0.7874;     % m
B = 0.4445;     % m
C = 0.0838;     % m
AA = 0.3454;    % m
BB = 0.1651;    % m

% Derived body dimensions
l_T = AA-C;             % m
l_N = B-AA;             % m
l_H = (A-C-l_T-l_N)/2;  % m

% Mass properties of torso and head
%m_T = 12.02+13.25;                     % kg (+/- 0.28) - torso (U+L)
%m_H = 3.73+0.91;                       % kg (+/- 0.14) - head and neck
m_T = 12.02+13.25+2*(1.18+0.90+0.28);   % kg (+/- 0.28) - torso (U+L), arms, hands
m_H = 3.73;                             % kg (+/- 0.05) - head
I_T = 0.35;                             % kg*m^2 - needs justification
I_H = (0.01^2)*(81.7*m_H-128.38);       % kg*m^2 - from Connor et al.

% Neck rotary spring/damper constants
% From Hoover & Meguid, assuming springs and dampers in series:
% UJSND: k = 5.39, c = 0
% JSND: k = 45.75, c = 0
% JDNS: k = 0, c = 2.54
% JSJD: k = 31.14, c = 10.29
k_N = 31.14;      % N*m/rad
c_N = 3.0;        % N*m*s/rad

% Airbag linear spring/damper constants
k_A = 500.0;      % N/m
c_A = 50.0;       % N*s/m

% Seat/seatbelt linear spring/damper constants
k_S = 5.0e4;      % N/m
c_S = 500.0;      % N*s/m

% Gravity
g = 9.81;         % m/s^2

%% Determine motion at the waist

% % NOTE: The commented-out stuff below can be used for notional prescribed
% % accelerations.
% 
% % Prescribed base acceleration at waist (notional)
% A_crash = 27*g;    % m/s^2 (peak acceleration)
% T_crash = 0.12;    % s (pulse duration)
% tspan = [0 0.3];   % s (total simulation time)
% x_W_ddot = @(t) -(t >= 0 & t <= T_crash).*A_crash.*sin(pi*t/T_crash).^2;
% y_W_ddot = @(t) zeros(size(t));
% 
% % Waist motion initial conditions (actual)
% x_W_0 = 0;              % m
% y_W_0 = 0;              % m
% x_W_dot_0 = 15.64;      % m/s (34.98 mph)
% y_W_dot_0 = 0;          % m/s
% 
% % Time vector
% t = linspace(tspan(1), tspan(2), 1000);
% 
% % Horizontal acceleration, velocity, and position
% x_W_ddot_vals = x_W_ddot(t);
% x_W_dot_vals = x_W_dot_0 + cumtrapz(t, x_W_ddot_vals);
% x_W_vals = x_W_0 + cumtrapz(t, x_W_dot_vals);
% 
% % Vertical acceleration, velocity, and position
% y_W_ddot_vals = y_W_ddot(t);
% y_W_dot_vals = y_W_dot_0 + cumtrapz(t, y_W_ddot_vals);
% y_W_vals = y_W_0 + cumtrapz(t, y_W_dot_vals);
% 
% % Interpolation functions for ode45
% x_W = @(tq) interp1(t, x_W_vals, tq, 'linear', 'extrap');
% x_W_dot = @(tq) interp1(t, x_W_dot_vals, tq, 'linear', 'extrap');
% y_W = @(tq) interp1(t, y_W_vals, tq, 'linear', 'extrap');
% y_W_dot = @(tq) interp1(t, y_W_dot_vals, tq, 'linear', 'extrap');

% Prescribed horizontal base acceleration at waist
% NOTE: The vertical component of acceleration (vehicle_az_proc) is quite
% small.
load('vehicle_accel_data.mat', ...
    'tdata_proc', 'vehicle_ax_proc', 'vehicle_az_proc');
t_data = tdata_proc;
ax_data = vehicle_ax_proc;
ay_data = vehicle_az_proc;

% Simulation time
tspan = [0 0.2999];

% Interpolation function for ode45
x_W_ddot = @(tq) interp1(t_data, ax_data, tq, 'linear', 0);

% No vertical acceleration
y_W_ddot = @(tq) interp1(t_data, ay_data, tq, 'linear', 0);

% Waist motion initial conditions
x_W_0 = 0;
y_W_0 = 0;
x_W_dot_0 = 15.64;
y_W_dot_0 = 0;

% Time vector
t = linspace(tspan(1), tspan(2), 1000);

% Horizontal acceleration, velocity, and position
x_W_ddot_vals = x_W_ddot(t);
x_W_dot_vals = x_W_dot_0 + cumtrapz(t, x_W_ddot_vals);
x_W_vals = x_W_0 + cumtrapz(t, x_W_dot_vals);

% Vertical acceleration, velocity, and position
y_W_ddot_vals = y_W_ddot(t);
y_W_dot_vals = y_W_dot_0 + cumtrapz(t, y_W_ddot_vals);
y_W_vals = y_W_0 + cumtrapz(t, y_W_dot_vals);

% Interpolation functions for ode45
x_W = @(tq) interp1(t, x_W_vals, tq, 'linear', 'extrap');
x_W_dot = @(tq) interp1(t, x_W_dot_vals, tq, 'linear', 'extrap');
y_W = @(tq) interp1(t, y_W_vals, tq, 'linear', 'extrap');
y_W_dot = @(tq) interp1(t, y_W_dot_vals, tq, 'linear', 'extrap');

%% Integrate EOMs

% Initial conditions
theta_0 = (90+1.7)*pi/180;
phi_0 = pi/2;
theta_dot_0 = 0;
phi_dot_0 = 0;
q0 = [theta_0;
      phi_0;
      theta_dot_0;
      phi_dot_0];

% Solve IVP
opts = odeset('MaxStep', 0.1, 'RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, q_sol] = ode45( ...
    @(t, q) EOMs_Lagrange( ...
        t, q, ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A, c_A, k_S, c_S, g, ...
        theta_0, phi_0, x_W_ddot(t), y_W_ddot(t)), ...
    tspan, q0);
    %tspan, q0, opts);

% Extract angles and angular velocities
theta = q_sol(:, 1);
phi = q_sol(:, 2);
theta_dot = q_sol(:, 3);
phi_dot = q_sol(:, 4);

% Trajectories of points of interest

% Waist
x_W_sol = x_W(t_sol);
y_W_sol = y_W(t_sol);

% Torso COM
x_T = x_W_sol + l_T*cos(theta);
y_T = y_W_sol + l_T*sin(theta);

% Base of neck
x_N = x_W_sol + (l_T + l_N)*cos(theta);
y_N = y_W_sol + (l_T + l_N)*sin(theta);

% Head COM
x_H = x_W_sol + (l_T + l_N)*cos(theta) + l_H*cos(phi);
y_H = y_W_sol + (l_T + l_N)*sin(theta) + l_H*sin(phi);

%% Head acceleration

% Evaluate angular accelerations from the equations of motion
theta_ddot = zeros(size(t_sol));
phi_ddot = zeros(size(t_sol));

for k = 1:length(t_sol)

    dqdt = EOMs_Lagrange( ...
        t_sol(k), q_sol(k,:)', ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A, c_A, k_S, c_S, g, ...
        theta_0, phi_0, ...
        x_W_ddot(t_sol(k)), y_W_ddot(t_sol(k)));

    theta_ddot(k) = dqdt(3);
    phi_ddot(k) = dqdt(4);

end

% Head acceleration from kinematics
x_H_ddot = x_W_ddot(t_sol) ...
    - (l_T + l_N)*cos(theta).*theta_dot.^2 ...
    - (l_T + l_N)*sin(theta).*theta_ddot ...
    - l_H*cos(phi).*phi_dot.^2 ...
    - l_H*sin(phi).*phi_ddot;

y_H_ddot = y_W_ddot(t_sol) ...
    - (l_T + l_N)*sin(theta).*theta_dot.^2 ...
    + (l_T + l_N)*cos(theta).*theta_ddot ...
    - l_H*sin(phi).*phi_dot.^2 ...
    + l_H*cos(phi).*phi_ddot;

%% Test data

load('head_accel_test_data.mat', 't001', 'a001', 't002', 'a002');

figure;
hold on;
plot(t001*1000, a001, 'LineWidth', 1.5, 'DisplayName', 'Test x Acceleration');
plot(t002*1000, a002, 'LineWidth', 1.5, 'DisplayName', 'Test y Acceleration');
xlabel('Time (ms)'); 
ylabel('Acceleration (g)'); 
title('Raw Test Accelerometer Data');
legend('Location', 'best'); grid on;

%% Simulation and test data comparison

a001_ms2 = a001*9.81;
a002_ms2 = a002*9.81;

figure;
hold on;
plot(t_sol*1000, x_H_ddot, 'LineWidth', 2, 'DisplayName', 'Simulated $\ddot{x}_H$');
plot(t_sol*1000, y_H_ddot, 'LineWidth', 2, 'DisplayName', 'Simulated $\ddot{y}_H$');
plot(t001*1000, a001_ms2, 'LineWidth', 1.5, 'DisplayName', 'Test x Acceleration');
plot(t002*1000, a002_ms2, 'LineWidth', 1.5, 'DisplayName', 'Test y Acceleration');
xlabel('Time (ms)'); 
ylabel('Acceleration (m/s^2)');
title('Simulated and Experimental Head Acceleration');
legend('Location', 'best', 'Interpreter', 'latex'); grid on;

%% Head acceleration

figure;
hold on;
plot(t_sol*1000, x_H_ddot, 'LineWidth', 2, 'DisplayName', '$\ddot{x}_H$');
plot(t_sol*1000, y_H_ddot, 'LineWidth', 2, 'DisplayName', '$\ddot{y}_H$');
xlabel('Time (ms)'); 
ylabel('Acceleration (m/s^2)');
title('Simulated Head Acceleration');
legend('Location', 'best', 'Interpreter', 'latex'); grid on;

%% Plotting results

% Prescribed waist acceleration
figure;
hold on;
plot(t*1000, x_W_ddot_vals, 'LineWidth', 2, 'DisplayName', '$\ddot{x}_W$', 'Color', 'k');
plot(t*1000, y_W_ddot_vals, 'LineWidth', 2, 'DisplayName', '$\ddot{y}_W$', 'Color', 'g');
xlabel('Time (ms)');
ylabel('Acceleration (m/s^2)')
title('Prescribed Waist Acceleration');
legend('Location', 'best', 'Interpreter', 'latex');
grid on;

% Angles of torso and head
figure;
hold on;
plot(t_sol*1000, theta*180/pi, 'LineWidth', 2, 'DisplayName', '\theta', 'Color', 'b');
plot(t_sol*1000, phi*180/pi, 'LineWidth', 2, 'DisplayName', '\phi', 'Color', 'r');
xlabel('Time (ms)');
ylabel('Angle (deg)');
legend('Location', 'best');
title('Torso and Head Angular Response');
grid on;

% X and Y of torso and head over time
figure;
hold on;
plot(t_sol*1000, x_T, 'LineWidth', 2, 'DisplayName', 'x_T', 'Color', 'b');
plot(t_sol*1000, y_T, 'LineWidth', 2, 'DisplayName', 'y_T', 'Color', 'b', 'LineStyle', '--');
plot(t_sol*1000, x_H, 'LineWidth', 2, 'DisplayName', 'x_H', 'Color', 'r');
plot(t_sol*1000, y_H, 'LineWidth', 2, 'DisplayName', 'y_H', 'Color', 'r', 'LineStyle', '--');
xlabel('Time (ms)');
ylabel('Displacement (m)');
legend('Location', 'best');
title('Torso and Head Angular Response');
grid on;

% Trajectories of torso and head
figure;
hold on;
plot(x_T, y_T, 'LineWidth', 2, 'DisplayName', 'Torso', 'Color', 'b');
plot(x_H, y_H, 'LineWidth', 2, 'DisplayName', 'Head', 'Color', 'r');
xlabel('x (m)');
ylabel('y (m)');
legend('Location', 'best');
title('Trajectory of Torso and Head');
axis equal;
grid on;

%% Animation of crash event
animate_motion(t_sol, x_W(t_sol), y_W(t_sol), ...
    x_T, y_T, x_N, y_N, x_H, y_H, ...
    theta, phi, l_T, l_N, l_H);

%% Parameter Study

% Parameter ranges
k_A_values = linspace(100, 1000, 10);
c_A_values = linspace(10, 150, 10);
k_S_values = linspace(1e4, 1e5, 10);
c_S_values = linspace(100, 1500, 10);

% Preallocate response metrics
peak_accel_kA = zeros(size(k_A_values));
peak_accel_cA = zeros(size(c_A_values));
peak_accel_kS = zeros(size(k_S_values));
peak_accel_cS = zeros(size(c_S_values));

peak_angle_kA = zeros(size(k_A_values));
peak_angle_cA = zeros(size(c_A_values));
peak_angle_kS = zeros(size(k_S_values));
peak_angle_cS = zeros(size(c_S_values));

%% Vary airbag stiffness

for i = 1:length(k_A_values)

    [t_p, q_p, xH_ddot_p, yH_ddot_p] = parameter_study( ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A_values(i), c_A, ...
        k_S, c_S, g, theta_0, phi_0, q0, ...
        tspan, x_W_ddot, y_W_ddot);

    a_H = sqrt(xH_ddot_p.^2 + yH_ddot_p.^2);

    peak_accel_kA(i) = max(a_H);
    peak_angle_kA(i) = max(abs(q_p(:,2) - phi_0))*180/pi;

end

%% Vary airbag damping

for i = 1:length(c_A_values)

    [t_p, q_p, xH_ddot_p, yH_ddot_p] = parameter_study( ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A, c_A_values(i), ...
        k_S, c_S, g, theta_0, phi_0, q0, ...
        tspan, x_W_ddot, y_W_ddot);

    a_H = sqrt(xH_ddot_p.^2 + yH_ddot_p.^2);

    peak_accel_cA(i) = max(a_H);
    peak_angle_cA(i) = max(abs(q_p(:,2) - phi_0))*180/pi;

end

%% Vary seatbelt stiffness

for i = 1:length(k_S_values)

    [t_p, q_p, xH_ddot_p, yH_ddot_p] = parameter_study( ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A, c_A, ...
        k_S_values(i), c_S, g, theta_0, phi_0, q0, ...
        tspan, x_W_ddot, y_W_ddot);

    a_H = sqrt(xH_ddot_p.^2 + yH_ddot_p.^2);

    peak_accel_kS(i) = max(a_H);
    peak_angle_kS(i) = max(abs(q_p(:,2) - phi_0))*180/pi;

end

%% Vary seatbelt damping

for i = 1:length(c_S_values)

    [t_p, q_p, xH_ddot_p, yH_ddot_p] = parameter_study( ...
        m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
        k_N, c_N, k_A, c_A, ...
        k_S, c_S_values(i), g, theta_0, phi_0, q0, ...
        tspan, x_W_ddot, y_W_ddot);

    a_H = sqrt(xH_ddot_p.^2 + yH_ddot_p.^2);

    peak_accel_cS(i) = max(a_H);
    peak_angle_cS(i) = max(abs(q_p(:,2) - phi_0))*180/pi;

end

%% Parameter Study Results

figure;

subplot(2,2,1);
plot(k_A_values, peak_accel_kA, 'o-', 'LineWidth', 2);
xlabel('$k_A$ (N/m)', 'Interpreter', 'latex');
ylabel('Peak Head Acceleration (m/s$^2$)', 'Interpreter', 'latex');
title('Effect of Airbag Stiffness');
grid on;

subplot(2,2,2);
plot(c_A_values, peak_accel_cA, 'o-', 'LineWidth', 2);
xlabel('$c_A$ (N$\cdot$s/m)', 'Interpreter', 'latex');
ylabel('Peak Head Acceleration (m/s$^2$)', 'Interpreter', 'latex');
title('Effect of Airbag Damping');
grid on;

subplot(2,2,3);
plot(k_S_values, peak_accel_kS, 'o-', 'LineWidth', 2);
xlabel('$k_S$ (N/m)', 'Interpreter', 'latex');
ylabel('Peak Head Acceleration (m/s$^2$)', 'Interpreter', 'latex');
title('Effect of Seatbelt Stiffness');
grid on;

subplot(2,2,4);
plot(c_S_values, peak_accel_cS, 'o-', 'LineWidth', 2);
xlabel('$c_S$ (N$\cdot$s/m)', 'Interpreter', 'latex');
ylabel('Peak Head Acceleration (m/s$^2$)', 'Interpreter', 'latex');
title('Effect of Seatbelt Damping');
grid on;

sgtitle('Parameter Study: Peak Head Acceleration');


%% FUNCTIONS

function dqdt = EOMs_Lagrange(~, q, m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
    k_N, c_N, k_A, c_A, k_S, c_S, g, theta_0, phi_0, x_W_ddot, y_W_ddot)

    % State variables
    theta = q(1);
    phi = q(2);
    theta_dot = q(3);
    phi_dot = q(4);

    % TODO: Reformat EOMs from ugly one-liners to less ugly multi-liners

    % theta double dot
    theta_ddot = ((2*theta_dot)*(l_H^2*m_H*c_A*(l_N + l_T)^2*cos(phi)^2 - l_H^2*((c_A - c_S)*l_T^2 + (2*l_N)*l_T*c_A + l_N^2*c_A)*m_H/2 + ((c_A + c_S)*l_T^2 + (2*l_N)*l_T*c_A + l_N^2*c_A)*I_H/2)*sin(theta)^2 + (-l_H^3*m_H*k_A*(l_N + l_T)*cos(phi)^3 + (2*m_H)*((m_H*theta_dot^2 - k_A/2)*(l_N + l_T)*cos(theta) + l_H*phi_dot*sin(phi)*c_A + k_A*(l_N + l_T)*cos(theta_0)/2 + l_H*cos(phi_0)*k_A/2 - m_H*x_W_ddot/2)*l_H^2*(l_N + l_T)*cos(phi)^2 - ((2*l_H)*m_H*theta_dot*sin(phi)*c_A*(l_N + l_T)*cos(theta) + l_H*m_H^2*(g + y_W_ddot)*sin(phi) - l_H^2*m_H^2*phi_dot^2 - m_H*phi_dot^2*I_H + I_H*k_A)*l_H*(l_N + l_T)*cos(phi) + (-l_H^2*theta_dot^2*(l_N + l_T)^2*m_H^2 - k_S*l_H^2*l_T^2*m_H - ((k_A + k_S)*l_T^2 + (2*k_A)*l_N*l_T + k_A*l_N^2)*I_H)*cos(theta) - l_H*(l_N + l_T)*((l_H^2*phi_dot*c_A - theta_dot*c_N + phi_dot*c_N - k_N*(theta - phi))*m_H - phi_dot*I_H*c_A)*sin(phi) + (k_S*l_H^2*l_T^2*m_H + ((k_A + k_S)*l_T^2 + (2*k_A)*l_N*l_T + k_A*l_N^2)*I_H)*cos(theta_0) + l_H*I_H*k_A*(l_N + l_T)*cos(phi_0) - ((l_H^2*l_T*m_T + I_H*(l_N + l_T))*m_H + I_H*l_T*m_T)*x_W_ddot)*sin(theta) + l_H^3*m_H*sin(phi)*cos(theta)*k_A*(l_N + l_T)*cos(phi)^2 - (2*m_H)*(sin(phi)*l_H*(m_H*theta_dot^2 - k_A/2)*(l_N + l_T)*cos(theta)^2 + (l_H^2*phi_dot*sin(phi)^2*c_A + (k_A*(l_N + l_T)*cos(theta_0) + l_H*cos(phi_0)*k_A - m_H*x_W_ddot)*l_H*sin(phi)/2 - theta_dot*c_N/2 + phi_dot*c_N/2 - k_N*(theta - phi)/2)*cos(theta) - l_H*m_H*theta_dot^2*sin(phi)*(l_N + l_T)/2)*l_H*(l_N + l_T)*cos(phi) + (l_H^2*m_H^2*(l_N + l_T)*(g + y_W_ddot)*sin(phi)^2 - l_H*m_H*phi_dot^2*(l_H^2*m_H + I_H)*(l_N + l_T)*sin(phi) + (g + y_W_ddot)*((l_H^2*l_T*m_T + I_H*(l_N + l_T))*m_H + I_H*l_T*m_T))*cos(theta) + (l_H^2*m_H + I_H)*(theta_dot*c_N - phi_dot*c_N + k_N*(theta - phi)))*1/((2*l_H^2)*m_H^2*cos(phi)*sin(phi)*sin(theta)*(l_N + l_T)^2*cos(theta) + (2*m_H^2)*(cos(theta)^2 - 1/2)*l_H^2*(l_N + l_T)^2*cos(phi)^2 - l_H^2*m_H^2*(l_N + l_T)^2*cos(theta)^2 + ((-l_T^2*m_T - I_T)*l_H^2 - I_H*(l_N + l_T)^2)*m_H - I_H*(l_T^2*m_T + I_T));

    % phi double dot
    phi_ddot = (-l_H*m_H*sin(phi)*k_A*(l_N + l_T)^3*cos(theta)^3 - m_H*(-l_H*phi_dot*c_A*(l_N + l_T)^2*sin(phi)^2 + (-theta_dot*c_A*(l_N + l_T)^3*sin(theta) - 2*(m_H*phi_dot^2 - k_A/2)*l_H*(l_N + l_T)^2*cos(phi) - k_A*(l_N + l_T)^3*cos(theta_0) - l_H*k_A*(l_N + l_T)^2*cos(phi_0) + ((l_N + l_T)^2*m_H + l_T^2*m_T)*x_W_ddot)*sin(phi) + (-((k_A + k_S)*l_T^2 + (2*k_A)*l_N*l_T + k_A*l_N^2)*(l_N + l_T)*sin(theta) + l_N*l_T*m_T*(g + y_W_ddot))*cos(phi))*l_H*cos(theta)^2 - l_H*(l_N + l_T)*((-k_S*l_T^2*m_H*sin(theta)^2 + m_H*(l_H*phi_dot*c_A*(l_N + l_T)*cos(phi) + (g + y_W_ddot)*((l_N + l_T)*m_H + l_T*m_T))*sin(theta) - theta_dot^2*(l_N + l_T)^2*m_H^2 - theta_dot^2*(l_T^2*m_T + I_T)*m_H + k_A*(l_T^2*m_T + I_T))*sin(phi) + m_H*(theta_dot*cos(phi)*((c_A + c_S)*l_T^2 + (2*l_N)*l_T*c_A + l_N^2*c_A)*sin(theta)^2 + (2*(m_H*phi_dot^2 - k_A/2)*l_H*(l_N + l_T)*cos(phi)^2 + (((k_A + k_S)*l_T^2 + (2*k_A)*l_N*l_T + k_A*l_N^2)*cos(theta_0) + cos(phi_0)*k_A*l_H*(l_N + l_T) - x_W_ddot*((l_N + l_T)*m_H + l_T*m_T))*cos(phi) - l_H*m_H*phi_dot^2*(l_N + l_T))*sin(theta) + cos(phi)*(theta_dot*c_N - phi_dot*c_N + k_N*(theta - phi))))*cos(theta) + phi_dot*l_H^2*((l_N + l_T)^2*m_H + (2*l_T^2)*m_T + 2*I_T)*c_A*sin(phi)^2 - l_H*(l_T^2*m_H*theta_dot*c_S*(l_N + l_T)*sin(theta)^3 + m_H*(k_S*l_T*(l_N + l_T)*cos(theta_0) - l_N*m_T*x_W_ddot)*l_T*sin(theta)^2 - ((l_T^2*theta_dot*c_A + (2*l_T)*theta_dot*c_A*l_N + theta_dot*c_A*l_N^2 - theta_dot*c_N + phi_dot*c_N - k_N*(theta - phi))*m_H + (2*theta_dot)*c_A*(l_T^2*m_T + I_T))*(l_N + l_T)*sin(theta) + l_H*(phi_dot^2*(l_N + l_T)^2*m_H^2 + k_A*(l_T^2*m_T + I_T))*cos(phi) - k_A*(l_T^2*m_T + I_T)*(l_N + l_T)*cos(theta_0) - l_H*k_A*(l_T^2*m_T + I_T)*cos(phi_0) + m_H*x_W_ddot*I_T)*sin(phi) + m_H*(g + y_W_ddot)*((l_N + l_T)^2*m_H + l_T^2*m_T)*l_H*cos(phi)*sin(theta)^2 - m_H*((l_N + l_T)^2*m_H + l_T^2*m_T + I_T)*theta_dot^2*l_H*cos(phi)*(l_N + l_T)*sin(theta) + l_H*m_H*I_T*(g + y_W_ddot)*cos(phi) - ((l_N + l_T)^2*m_H + l_T^2*m_T + I_T)*(theta_dot*c_N - phi_dot*c_N + k_N*(theta - phi)))*1/(l_H^2*m_H^2*cos(phi)^2*(l_N + l_T)^2*cos(theta)^2 + (2*l_H^2)*m_H^2*cos(phi)*sin(phi)*sin(theta)*(l_N + l_T)^2*cos(theta) + l_H^2*m_H^2*sin(theta)^2*(l_N + l_T)^2*sin(phi)^2 - ((l_N + l_T)^2*m_H + l_T^2*m_T + I_T)*(l_H^2*m_H + I_H));

    % State derivatives
    dqdt = [
        theta_dot
        phi_dot
        theta_ddot
        phi_ddot
    ];

end


function animate_motion(t_sol, x_W, y_W, x_T, y_T, x_N, y_N, x_H, y_H, ...
    theta, phi, l_T, l_N, l_H)
    
    % Dimensions of ellipses to represent body parts
    % NOTE: These dimensions are purely visual
    torso_length = 0.30;   % m
    torso_width = 0.25;    % m
    head_length = 0.20;    % m
    head_width = 0.15;     % m
    
    % Initialize figure
    figure;
    hold on;
    grid on;
    axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    title('Crash Event Animation');
    
    % Axis limits
    x_margin = 0.5;
    y_margin = 0.5;
    xlim([min(x_W)-x_margin, max(x_W)+l_T+l_N+l_H+x_margin]);
    ylim([min(y_W)-y_margin, max(y_W)+l_T+l_N+l_H+y_margin]);
    
    % Create body segments
    h_torso = plot([0 0], [0 0], 'b-', 'LineWidth', 3);
    h_neck = plot([0 0], [0 0], 'k-', 'LineWidth', 3);
    h_head = plot([0 0], [0 0], 'r-', 'LineWidth', 3);
    
    % Create joints and COM points
    h_waist = plot(0, 0, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    h_T = plot(0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
    h_N = plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    h_H = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    
    % Create transparent torso ellipse
    [ex, ey] = ellipse_points(torso_length/2, torso_width/2, 0, 0, 0);
    h_torso_body = patch(ex, ey, 'b', 'FaceAlpha', 0.25, ...
        'EdgeColor', 'b', 'LineWidth', 1.5);
    
    % Create transparent head ellipse
    [ex, ey] = ellipse_points(head_length/2, head_width/2, 0, 0, 0);
    h_head_body = patch(ex, ey, 'r', 'FaceAlpha', 0.25, ...
        'EdgeColor', 'r', 'LineWidth', 1.5);
    
    % Animation
    for k = 1:length(t_sol)
    
        % Update torso, neck, and head segments
        set(h_torso, 'XData', [x_W(k) x_T(k)], 'YData', [y_W(k) y_T(k)]);
        set(h_neck, 'XData', [x_T(k) x_N(k)], 'YData', [y_T(k) y_N(k)]);
        set(h_head, 'XData', [x_N(k) x_H(k)], 'YData', [y_N(k) y_H(k)]);
    
        % Update joints
        set(h_waist, 'XData', x_W(k), 'YData', y_W(k));
        set(h_T, 'XData', x_T(k), 'YData', y_T(k));
        set(h_N, 'XData', x_N(k), 'YData', y_N(k));
        set(h_H, 'XData', x_H(k), 'YData', y_H(k));
    
        % Update torso ellipse
        [ex, ey] = ellipse_points(torso_length/2, torso_width/2, ...
            x_T(k), y_T(k), theta(k));
        set(h_torso_body, 'XData', ex, 'YData', ey);
    
        % Update head ellipse
        [ex, ey] = ellipse_points(head_length/2, head_width/2, ...
            x_H(k), y_H(k), phi(k));
        set(h_head_body, 'XData', ex, 'YData', ey);
    
        % Time display
        title(sprintf('Crash Event (t = %.3f s)', t_sol(k)));
        drawnow;
    
        pause(0.01);
    end

end


function [x, y] = ellipse_points(a, b, x0, y0, angle)

    t = linspace(0, 2*pi, 50);

    % Ellipse before rotation
    x_local = a*cos(t);
    y_local = b*sin(t);

    % Rotation matrix
    x = x0 + x_local*cos(angle) - y_local*sin(angle);
    y = y0 + x_local*sin(angle) + y_local*cos(angle);

end

function [t_sol, q_sol, x_H_ddot, y_H_ddot] = parameter_study( ...
    m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
    k_N, c_N, k_A, c_A, k_S, c_S, ...
    g, theta_0, phi_0, q0, ...
    tspan, x_W_ddot, y_W_ddot)

    % Integrate equations of motion
    [t_sol, q_sol] = ode45( ...
        @(t,q) EOMs_Lagrange( ...
            t, q, ...
            m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
            k_N, c_N, k_A, c_A, k_S, c_S, ...
            g, theta_0, phi_0, ...
            x_W_ddot(t), y_W_ddot(t)), ...
        tspan, q0);

    % Extract states
    theta = q_sol(:,1);
    phi = q_sol(:,2);
    theta_dot = q_sol(:,3);
    phi_dot = q_sol(:,4);

    % Evaluate angular accelerations
    theta_ddot = zeros(size(t_sol));
    phi_ddot = zeros(size(t_sol));

    for k = 1:length(t_sol)

        dqdt = EOMs_Lagrange( ...
            t_sol(k), q_sol(k,:)', ...
            m_T, m_H, I_T, I_H, l_T, l_N, l_H, ...
            k_N, c_N, k_A, c_A, k_S, c_S, ...
            g, theta_0, phi_0, ...
            x_W_ddot(t_sol(k)), y_W_ddot(t_sol(k)));

        theta_ddot(k) = dqdt(3);
        phi_ddot(k) = dqdt(4);

    end

    % Head acceleration
    x_H_ddot = x_W_ddot(t_sol) ...
        - (l_T + l_N)*cos(theta).*theta_dot.^2 ...
        - (l_T + l_N)*sin(theta).*theta_ddot ...
        - l_H*cos(phi).*phi_dot.^2 ...
        - l_H*sin(phi).*phi_ddot;

    y_H_ddot = y_W_ddot(t_sol) ...
        - (l_T + l_N)*sin(theta).*theta_dot.^2 ...
        + (l_T + l_N)*cos(theta).*theta_ddot ...
        - l_H*sin(phi).*phi_dot.^2 ...
        + l_H*cos(phi).*phi_ddot;

end