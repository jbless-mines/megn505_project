% MEGN505
% Crash test data analysis

%% Set-up

clc; clear all; close all;

% Plot settings
set(groot, 'defaultAxesFontSize', 14)
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.4)
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.4)

% Import crash test datafile
load('subaru_legacy_passenger_data.mat');

% Checking if the angular displacement is primarily positive or negative
beta = cumtrapz(tdata,omega2);
figure()
clf;
plot(tdata*1000,beta,'k-','LineWidth',2);
xlabel('Time (ms)')
ylabel('Estimated head angular displacement (deg)')
grid()

% If beta is primarily positive, switch signs of omega values
%omega1 = -omega1;
%omega2 = -omega2;
%omega3 = -omega3;

% Converting acceleration data from g's to m/s^2
g = 9.81;   % m/s^2
a1 = a1*g;
a2 = a2*g;
a3 = a3*g;

%% Calculating Euler angles using ODE45

v0 = [0;0;0];       % ICs
tol = 1e-8;         % simulation tolerance value
opts = odeset('abstol',tol,'reltol',tol);  % setting tolerances
[t,v] = ode45(@(t,v)orientation_from_gyroscope(t,v,omega1,omega2,omega3,tdata),tdata,v0,opts);

psi = v(:,1);
theta = v(:,2);
phi = v(:,3);

figure()
clf;
subplot(3,1,1)
plot(tdata*1000,psi,'k-','LineWidth',2)
ylim([-100 100])
ylabel('\psi (deg)')
xlabel('t (ms)')
grid()
subplot(3,1,2)
plot(tdata*1000,theta,'k-','LineWidth',2)
ylim([-100 100])
ylabel('\theta (deg)')
xlabel('t (ms)')
grid()
subplot(3,1,3)
plot(tdata*1000,phi,'k-','LineWidth',2)
ylim([-100 100])
ylabel('\phi (deg)')
xlabel('t (ms)')
grid()
sgtitle('Head Orientation','FontSize',16,'FontWeight','bold')

%% Transforming acceleration data

for i=1:size(tdata,1)
    a1_now = a1(i);
    a2_now = a2(i);
    a3_now = a3(i);
    psi_now = psi(i);
    theta_now = theta(i);
    phi_now = phi(i);
    R_b_I = [cosd(psi_now)*cosd(theta_now),cosd(psi_now)*sind(theta_now)*sind(phi_now)-sind(psi_now)*cosd(phi_now),cosd(psi_now)*sind(theta_now)*cosd(phi_now)+sind(psi_now)*sind(phi_now);
             sind(psi_now)*cosd(theta_now),sind(psi_now)*sind(theta_now)*sind(phi_now)+cosd(psi_now)*cosd(phi_now),sind(psi_now)*sind(theta_now)*cosd(phi_now)-cosd(psi_now)*sind(phi_now);
             -sind(theta_now),cosd(theta_now)*sind(phi_now),cosd(theta_now)*cosd(phi_now)];
    aI(:,i) = R_b_I*[a1_now;a2_now;a3_now];
end

aX = aI(1,:)';
aY = aI(2,:)';
aZ = aI(3,:)';

%% Calculating displacements using ODE45

v0 = [0;0;0;35/2.237;0;0];  % ICs (mph converted to m/s)
tol = 1e-8;         % simulation tolerance value
opts = odeset('abstol',tol,'reltol',tol);  % setting tolerances
[t,v] = ode45(@(t,v)disp_from_accel(t,v,aX,aY,aZ,tdata),tdata,v0,opts);

X = v(:,1);
Y = v(:,2);
Z = v(:,3);

figure()
clf;
subplot(3,1,1)
plot(tdata*1000,X,'k-','LineWidth',2)
ylim([0 2])
ylabel('X (m)')
xlabel('t (ms)')
grid()
subplot(3,1,2)
plot(tdata*1000,Y,'k-','LineWidth',2)
ylim([-.1 .2])
ylabel('Y (m)')
xlabel('t (ms)')
grid()
subplot(3,1,3)
plot(tdata*1000,Z,'k-','LineWidth',2)
ylim([-.2 .4])
ylabel('Z (m)')
xlabel('t (ms)')
grid()
sgtitle('Head Position','FontSize',16,'FontWeight','bold')

%% Animation of head motion

% NOTE: View options include 'default', 'passenger', and 'driver'
animate_head_motion(X,Y,Z,psi*pi/180,theta*pi/180,phi*pi/180,3,2,1,'default')

%% Import and post-process vehicle acceleration data

% Import data
% NOTE: ax and axr are intended to be redundant measurements. az is
% vertical, positive downward (hence the -1 multiplier applied later).
vehicle_ax = readmatrix('v10150.094', 'FileType', 'text', 'Delimiter', '\t');
vehicle_axr = readmatrix('v10150.100', 'FileType', 'text', 'Delimiter', '\t');
vehicle_az = readmatrix('v10150.098', 'FileType', 'text', 'Delimiter', '\t');
vehicle_ax = vehicle_ax(:, 2)*g;
vehicle_axr = vehicle_axr(:, 2)*g;
vehicle_az = -vehicle_az(:, 2)*g;

% Process data
vehicle_ax_mean = 0.5*(vehicle_ax+vehicle_axr);
window = 101;
shift = 15;
vehicle_ax_proc = movmean(vehicle_ax_mean, window);
vehicle_ax_proc = [zeros(shift, 1); vehicle_ax_proc(1:end-shift)];
vehicle_ax_proc(tdata < 0 | tdata > 0.1) = 0;

% Process data
vehicle_az_proc = movmean(vehicle_az, window);
vehicle_az_proc = [zeros(shift, 1); vehicle_az_proc(1:end-shift)];
vehicle_az_proc(tdata < 0 | tdata > 0.1) = 0;

% Keep only time >= 0
idx = tdata >= 0;
tdata_proc = tdata(idx);
vehicle_ax_proc = vehicle_ax_proc(idx);
vehicle_az_proc = vehicle_az_proc(idx);

% Plot data
figure();
hold on;
plot(tdata*1000, vehicle_ax_mean, 'r', 'LineWidth', 2, 'DisplayName', 'Original')
plot(tdata_proc*1000, vehicle_ax_proc, 'k', 'LineWidth', 2, 'DisplayName', 'Processed')
title('Vehicle Longitudinal Acceleration')
xlabel('Time (ms)')
ylabel('Acceleration (m/s^2)')
legend('Location', 'best')
grid()

% Plot data
figure();
hold on;
plot(tdata*1000, vehicle_az, 'r', 'LineWidth', 2, 'DisplayName', 'Original')
plot(tdata_proc*1000, vehicle_az_proc, 'k', 'LineWidth', 2, 'DisplayName', 'Processed')
title('Vehicle Vertical Acceleration')
xlabel('Time (ms)')
ylabel('Acceleration (m/s^2)')
legend('Location', 'best')
grid()

% Export
save('vehicle_accel_data.mat', 'tdata_proc', 'vehicle_ax_proc', 'vehicle_az_proc');
