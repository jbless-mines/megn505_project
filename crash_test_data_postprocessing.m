% MEGN505
% Final Project
% John Blessinger & Adrian Davis
% crash_test_test_data_postprocessing.m

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
animate_head_motion(X,Y,Z,psi*pi/180,theta*pi/180,phi*pi/180,3,2,1,tdata,'passenger')

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

%% Export processed data

start_idx = find(tdata == tdata_proc(1));
X_proc = X(start_idx:end);
Y_proc = Y(start_idx:end);
Z_proc = Z(start_idx:end);
psi_proc = psi(start_idx:end);
theta_proc = theta(start_idx:end);
phi_proc = phi(start_idx:end);

save('head_position_data.mat', ...
    'tdata_proc', ...
    'X_proc', 'Y_proc', 'Z_proc', ...
    'phi_proc', 'theta_proc', 'phi_proc')
save('vehicle_accel_data.mat', ...
    'tdata_proc', 'vehicle_ax_proc', 'vehicle_az_proc');

%% FUNCTIONS

function animate_head_motion(X, Y, Z, psi, theta, phi, r1, r2, r3, t, pov)

%  t should contain the time corresponding to each frame, in seconds.
%  For example: t = (0:length(X)-1)/1000;  % 1000 Hz data

%  For simplicity, represent the head as a sphere:

Rhead = 0.1;                              %  m
[xhead, yhead, zhead] = sphere(36);       %  m
xhead = Rhead*xhead;                       %  m
yhead = Rhead*yhead;                       %  m
zhead = Rhead*zhead;                       %  m

%  Include two small spheres for eyes:

Reye = 0.015;                              %  m
[xeye, yeye, zeye] = sphere(36);           %  m
xeye = Reye*xeye;                          %  m
yeye = Reye*yeye;                          %  m
zeye = Reye*zeye;                          %  m  

alt = 120*(pi/180);                        %  rad
azi = 30*(pi/180);                         %  rad

xrighteye = xeye + Rhead*sin(alt)*cos(+azi); 
yrighteye = yeye + Rhead*sin(alt)*sin(+azi); 
zrighteye = zeye + Rhead*cos(alt);                   

xlefteye = xeye + Rhead*sin(alt)*cos(+azi);         
ylefteye = yeye + Rhead*sin(alt)*sin(-azi);         
zlefteye = zeye + Rhead*cos(alt);                   

%  Create a circle for drawing lines of longitude and latitude on the head:

angle = linspace(0, 2*pi, 40)';             %  rad
circle1 = Rhead*cos(angle);                 %  m
circle2 = Rhead*sin(angle);                 %  m

circlexy = [circle1, circle2, 0*circle2];   %  m
circleyz = [0*circle1, circle1, circle2];   %  m
circlexz = [circle1, 0*circle1, circle2];   %  m

%  Set up the figure window:       

figure(101)
set(gcf, 'color', 'w')
plot3(X(1), Y(1), Z(1))
xlabel('\itX\rm (m)')
set(gca, 'xdir', 'reverse')
ylabel('\itY\rm (m)')
if strcmpi(pov, 'default')
    zlabel('\itZ\rm (m)            ', 'rotation', 0)
else
    zlabel('\itZ\rm (m)', 'rotation', 0)
end
set(gca, 'zdir', 'reverse')
axis equal
xlim([min(X)-1.2*Rhead, max(X)+1.2*Rhead])
ylim([min(Y)-1.2*Rhead, max(Y)+1.2*Rhead])
zlim([min(Z)-1.2*Rhead, max(Z)+1.2*Rhead])
grid on

%  Set the view based on the user input: 

if strcmp(pov, 'default') == 1
   view(-40, 20)
elseif strcmp(pov, 'driver') == 1
   view(0, 0)
elseif strcmp(pov, 'passenger') == 1
   view(180, 0)     
else
   view(-40, 20)
end

%  Add a title that will be updated during the animation:

time_title = title(sprintf('Time = %.4f s', t(1)));

%  Assemble the components to draw the head:

head = hgtransform;

dummy_head = surf('xdata', xhead, 'ydata', yhead, 'zdata', zhead, ...
                  'edgecolor', 'none', 'facecolor', [1, 0.85, 0.70], ...
                  'facealpha', 1, 'edgealpha', 1, 'parent', head);

right_eye = surf('xdata', xrighteye, 'ydata', yrighteye, ...
                 'zdata', zrighteye, 'edgecolor', 'none', ...
                 'facecolor', 'k', 'facealpha', 1, 'edgealpha', 1, ...
                 'parent', head);
             
left_eye = surf('xdata', xlefteye, 'ydata', ylefteye, ...
                'zdata', zlefteye, 'edgecolor', 'none', ...
                'facecolor', 'k', 'facealpha', 1, 'edgealpha', 1, ...
                'parent', head);

circle_xy = line('xdata', circlexy(:,1), 'ydata', circlexy(:,2), ...
                 'zdata', circlexy(:,3), 'color', 'b', 'linewidth', 3, ...
                 'parent', head);
             
circle_yz = line('xdata', circleyz(:,1), 'ydata', circleyz(:,2), ...
                 'zdata', circleyz(:,3), 'color', 'r', 'linewidth', 3, ...
                 'parent', head);
             
circle_xz = line('xdata', circlexz(:,1), 'ydata', circlexz(:,2), ...
                 'zdata', circlexz(:,3), 'color', 'k', 'linewidth', 3, ...
                 'parent', head);

%  Use the provided Euler angle sequence to specify the appropriate axes
%  of rotation:

sequence = [r1, r2, r3];

for k = 1:3
    if sequence(k) == 1
       rotation_axes(:,k) = [1, 0, 0]';
    elseif sequence(k) == 2
       rotation_axes(:,k) = [0, 1, 0]';
    else
       rotation_axes(:,k) = [0, 0, 1]';
    end
end

%  Display the head in its initial position and orientation:

head.Matrix = makehgtform('translate', [X(1), Y(1), Z(1)]', ...
                          'axisrotate', rotation_axes(:,1), psi(1), ...
                          'axisrotate', rotation_axes(:,2), theta(1), ...
                          'axisrotate', rotation_axes(:,3), phi(1));

drawnow

%  Animate the head's motion:

pause(2)

% animation = VideoWriter(strcat('head-motion-', pov, '-view.avi'));
% animation.FrameRate = 100;
% open(animation);

for k = 1:5:length(X)

    head.Matrix = makehgtform('translate', [X(k), Y(k), Z(k)]', ...
                              'axisrotate', rotation_axes(:,1), psi(k), ...
                              'axisrotate', rotation_axes(:,2), theta(k), ...
                              'axisrotate', rotation_axes(:,3), phi(k));

    % Update the time displayed in the title:
    time_title.String = sprintf('Time = %.4f s', t(k));

    drawnow

    % writeVideo(animation, getframe(gcf));

end

% close(animation);

end

function vdot = disp_from_accel(t,v,aX,aY,aZ,tdata)

aX_now = interp1(tdata,aX,t);
aY_now = interp1(tdata,aY,t);
aZ_now = interp1(tdata,aZ,t);

xdot = v(4);
ydot = v(5);
zdot = v(6);

vdot = [v(4);
        v(5);
        v(6);
        aX_now;
        aY_now;
        aZ_now];

end


function vdot = orientation_from_gyroscope(t,v,omega1,omega2,omega3,tdata)

om1_now = interp1(tdata,omega1,t);
om2_now = interp1(tdata,omega2,t);
om3_now = interp1(tdata,omega3,t);

psi = v(1);
theta = v(2);
phi = v(3);

vdot = [om2_now*sind(phi)/cosd(theta) + om3_now*cosd(phi)/cosd(theta);
        om2_now*cosd(phi) - om3_now*sind(phi);
        om1_now + om2_now*sind(theta)*sind(phi)/cosd(theta) + om3_now*sind(theta)*cosd(phi)/cosd(theta)];

end