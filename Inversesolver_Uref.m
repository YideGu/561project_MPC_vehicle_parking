function U_ref=Inversesolver_Uref(Y_ref_in,dY_ref, dt)  

%{
num=604
Y_ref_in=Y_ref(:,num);
dY_ref=Y_ref_in-Y_ref(:,num-1);
%}
parametersetup
syms d delta
 
delta_val=double (vpasolve(dY_ref(6)/dt==1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*Y_ref_in(6) + Y_ref_in(5),Y_ref_in(4))+delta))))*l_f*cos(delta)- (D_r*sin(C_r*atan(B_r*(atan2(l_r*Y_ref_in(6) - Y_ref_in(5),Y_ref_in(4))))))*l_r),delta));

eqn2=subs(dY_ref(4)/dt==1/m*((Cm1*d-Cm2*d*Y_ref_in(4)-Cr0-Cr2*Y_ref_in(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*Y_ref_in(6) + Y_ref_in(5),Y_ref_in(4))+delta))))*sin(delta) + m*Y_ref_in(5)*Y_ref_in(6)),delta,delta_val);
d_val=double(vpasolve(eqn2,d));
if isempty(delta_val)
U_ref=[0, 0];
else
U_ref=[delta_val, d_val]
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