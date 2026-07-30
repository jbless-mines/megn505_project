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