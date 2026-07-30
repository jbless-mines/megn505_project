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