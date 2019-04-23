function [xfinal,yfinal,phifinal, v_xfinal,v_yfinal,omegafinal] = dynamic_model(dt,IC,u)

%{
K=0.1;
K_n=10;
G=1.0;
G_n=115;

u=[[G*ones(1,G_n), zeros(1,201-G_n)];  0.001*ones(1,201) ];
T=6,IC=[8.87175730240378;10.2000000000000;0;3.2;0;0]
%}

parametersetup
n=length(u);
px=zeros(n,1);
py=zeros(n,1);
phi=zeros(n,1);
v_x=zeros(n,1);
v_y=zeros(n,1);
omega=zeros(n,1);

                
 
odefun=@(t,x, u) [
    x(4)*cos(x(3))-x(5)*sin(x(3));
    x(4)*sin(x(3))+x(5)*cos(x(3));
    x(6);
    1/m*((Cm1*u(2)-Cm2*u(2)*x(4)-Cr0-Cr2*x(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*sin(u(1)) + m*x(5)*x(6));
    1/m*(D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))) + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*cos(u(1)) - m*x(4)*x(6));
    1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*l_f*cos(u(1))- (D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))))*l_r)
    ];
  px(1)=IC(1);
 py(1)=IC(2);
  phi(1)=IC(3);
 v_x(1)=IC(4);
v_y(1)=IC(5);
omega(1)=IC(6);

for i =1:1:n
[temp_t, temp_x] = ode23(@(t,x) odefun(t, x, [u(1,i) u(2,i)]), [0 dt], [px(i) py(i) phi(i) v_x(i) v_y(i) omega(i)]);
px(i+1) = temp_x(end,1);
py(i+1) = temp_x(end,2);
phi(i+1) = temp_x(end,3);
v_x(i+1) = temp_x(end,4);
v_y(i+1) = temp_x(end,5);
omega(i+1) = temp_x(end,6);


%figure

end




xfinal= px(2:n+1)';
yfinal =  py(2:n+1)';
phifinal =  phi(2:n+1)';
v_xfinal= v_x(2:n+1)';
v_yfinal =  v_y(2:n+1)';
omegafinal= omega(2:n+1)';

