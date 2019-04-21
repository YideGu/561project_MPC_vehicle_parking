function d_val=Inversesolver_d(u,Y_ref_in, dt)  
%{
u=[3.1 0]
Y_ref_in=[6.71006690135136;10.2000000000000;0;3.00191139677196;0;0];
dt=0.01;
%}
parametersetup

d_val=((((u(1)-Y_ref_in(4))/dt*m + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*Y_ref_in(6) + Y_ref_in(5),Y_ref_in(4))+abs(u(2))))))*sin(abs(u(2))) - m*Y_ref_in(5)*Y_ref_in(6)))+Cr0+Cr2*Y_ref_in(4)^2)/(Cm1-Cm2*Y_ref_in(4));
end  
%{
clear all
i=1
for delta=-0.5:0.01:0.5
x(i)=delta
y(i)=(3435758140106421*sin(2*atan(13*delta - 82538129168911975/144115188075855872))*cos(delta))/789724226650112 + 7737333415762779/6317793813200896
;
i=i+1;
end
plot(x,y)
%}