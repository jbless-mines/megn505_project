function animate_head_motion(X, Y, Z, psi, theta, phi, r1, r2, r3, pov)

%  For simplicity, represent the head as a sphere:

Rhead = 0.1;                             	%  m
[xhead, yhead, zhead] = sphere(36);     	%  m
xhead = Rhead*xhead;                     	%  m
yhead = Rhead*yhead;                    	%  m
zhead = Rhead*zhead;                        %  m

%  Include two small spheres for eyes:

Reye = 0.015;                               %  m
[xeye, yeye, zeye] = sphere(36);            %  m
xeye = Reye*xeye;                           %  m
yeye = Reye*yeye;                           %  m
zeye = Reye*zeye;                           %  m  

alt = 120*(pi/180);         %  rad
azi = 30*(pi/180);          %  rad

xrighteye = xeye + Rhead*sin(alt)*cos(+azi);        %  m
yrighteye = yeye + Rhead*sin(alt)*sin(+azi);        %  m
zrighteye = zeye + Rhead*cos(alt);                  %  m

xlefteye = xeye + Rhead*sin(alt)*cos(+azi);         %  m
ylefteye = yeye + Rhead*sin(alt)*sin(-azi);         %  m
zlefteye = zeye + Rhead*cos(alt);                   %  m
       
%  Create a circle for drawing lines of longitude and latitude on the head
%  to better visualize the head's change in orientation:

angle = linspace(0, 2*pi, 40)';             %  rad
circle1 = Rhead*cos(angle);                 %  m
circle2 = Rhead*sin(angle);             	%  m

circlexy = [circle1, circle2, 0*circle2];       	%  m
circleyz = [0*circle1, circle1, circle2];           %  m
circlexz = [circle1, 0*circle1, circle2];           %  m

%  Set up the figure window:       

figure(101)
set(gcf, 'color', 'w')
plot3(X(1), Y(1), Z(1))
xlabel('\itX\rm (m)')
set(gca, 'xdir', 'reverse')
ylabel('\itY\rm (m)')
zlabel('\itZ\rm (m)            ', 'rotation', 0)
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

%  Use the provided Euler angle sequence to specify the appropriate axes of
%  rotation:

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

%  Animate the head's motion by updating the figure with its current
%  location and orientation:

pause(2)

% animation = VideoWriter(strcat('head-motion-', pov, '-view.avi'));
% animation.FrameRate = 100;
% open(animation);

for k = 1:5:length(X)
    head.Matrix = makehgtform('translate', [X(k), Y(k), Z(k)]', ...
                              'axisrotate', rotation_axes(:,1), psi(k), ...
                              'axisrotate', rotation_axes(:,2), theta(k), ...
                              'axisrotate', rotation_axes(:,3), phi(k));
    drawnow
    % writeVideo(animation, getframe(gcf));
end

% close(animation);