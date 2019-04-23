function [xfinal,yfinal,phifinal] = kinematic_model(T,IC,b,L,u)

n=length(u);

dt=T/n;





odefun = @(t, x, u) [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
                  u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
                  u(1)*tan(u(2))/L] ;
              px(1)=IC(1);
              py(1)=IC(2);
              phi(1)=IC(3);
for i =1:1:T/dt
[temp_t, temp_x] = ode45(@(t,x) odefun(t, x, [u(1,i) u(2,i)]), [0 dt], [px(i) py(i) phi(i)]);
px(i+1) = temp_x(end,1);
py(i+1) = temp_x(end,2);
phi(i+1) = temp_x(end,3);
plot(px,py);
%figure
%plot(x, y);
end
xfinal= px(2:n+1);
yfinal =  py(2:n+1);
phifinal =  phi(2:n+1);

              
          