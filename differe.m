function dx=differe(x,delta,d)
persistent dxfunc
if isempty(dxfunc)
parametersetup;
dxfunc =@(x,delta,d)[
    x(4)*cos(x(3))-x(5)*sin(x(3));
    x(4)*sin(x(3))+x(5)*cos(x(3));
    x(6);
    1/m*((Cm1*d-Cm2*d*x(4)-Cr0-Cr2*x(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*sin(delta) + m*x(5)*x(6));
    1/m*(D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))) + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*cos(delta) - m*x(4)*x(6));
    1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*l_f*cos(delta)- (D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))))*l_r)
    ];
end

dx=dxfunc(x,delta,d);

