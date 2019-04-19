function [x,y,phi] = kinematic_model(T,IC,b,L,u,dt)

odefun = @(t, x, u) [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
                  u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
                  u(1)*tan(u(2))/L] ;

[temp_t, temp_x] = ode45(@(t,x) odefun(t, x, u), [0 T], IC);
x = temp_x(:,1);
y = temp_x(:,2);
phi = temp_x(:,3);
%plot(temp_t,temp_x);
%figure
%plot(x, y);
x = interp1(temp_t,x,0:dt:T)';
y = interp1(temp_t,y,0:dt:T)';
phi = interp1(temp_t,phi,0:dt:T)';

              
          