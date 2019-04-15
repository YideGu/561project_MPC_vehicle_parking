function [t,statenew]=trajectory(state,u,pf,p,T)
%%%%%%%%%%%%%%% 
% state vector: x
% input vector: u which is fixed during T
% T, discrete time 
% p,pf, pamameter vector
%%
%vecotr parameter of force
Cm1=pf(1);
Cm2=pf(2);
Cr0=pf(3);
Cr2=pf(4);
B_r=pf(5);
C_r=pf(6);
D_r=pf(7);
B_f=pf(8);
C_f=pf(9);
D_f=pf(10);
%parameter of system
m=p(1);
Iz=p(2);
l_f=p(3);
l_r=p(4);
%input 
delta=u(1);
d=u(2)   ;
% state 


%F_fy = D_f*sin(C_f*atan(B_f*( -atan2(l_f*omega + vy,vx)+delta)));
%F_ry = D_r*sin(C_r*atan(B_r*(atan2(l_r*omega - vy,vx))));
%F_rx = (Cm1*d-Cm2*d*vx-Cr0-Cr2*vx^2);


% nonlinear function
df=@(t,x) [
    x(4)*cos(x(3))-x(5)*sin(x(3));
    x(4)*sin(x(3))+x(5)*cos(x(3));
    x(6);
    1/m*((Cm1*d-Cm2*d*x(4)-Cr0-Cr2*x(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*sin(delta) + m*x(5)*x(6));
    1/m*(D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))) + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*cos(delta) - m*x(4)*x(6));
    1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+delta))))*l_f*cos(delta)- (D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))))*l_r)
    ];

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,statenew]=ode45(df,[0 T],state, opts);
statenew= statenew + [0,0,0,normrnd(0,0.01),normrnd(0,0.01),normrnd(0,0.01)];
end