clear;
clc;
%% parameters
%wheelbases
% bounding box 4.5m x 1.7m for compact car

parametersetup;
L=2.7;
%distance from rear wheel to center of mass
b=1.35;
%open Kalman filter: =1; or not: else value
 Kalman_on=0;
%time discretization
dt=0.01;
%time span We make the simulation in 18 seconds
T=0:dt:18;
x_store=zeros(3,length(T));
x_store_d=zeros(6,length(T));
X0 = [12 10 pi]';
X0_d = [12 10 pi 2 0 0]';
Y_ref(:,1)=[12 10 pi]';
  x_store(:,1)=X0;
    x_store_d(:,1)=X0_d;
%% load reference trajectory


U_ref(1:2,1:600)=[ones(1,600) * 2; zeros(1,600)];

[x,y,phi] = kinematic_model(6,Y_ref(:,1),b,L,U_ref(1:2,1:600));

Y_ref(1:3,2:601)=[x;y;phi];

U_ref(1:2,601:1201)=[ones(1,601) * 1; ones(1,601) * 0.8];

[x,y,phi] = kinematic_model(6,Y_ref(:,601),b,L,U_ref(1:2,601:1201));

Y_ref(1:3,601:1201)=[x;y;phi];

U_ref(1:2,1201:1801)=[ones(1,601) * 2; ones(1,601) * 0]; %0.6

[x,y,phi] = kinematic_model(6,Y_ref(:,1201),b,L,U_ref(1:2,1201:1801));

Y_ref(1:3,1201:1801)=[x;y;phi];











%the reference trajectory Yref is given in a 3x601 double

%% 2.1 Discrete-time A and B matrices
%these are the system linearized in discrete time about the reference
%trajectory i.e. x(i+1)=A_i*x_i+B_i*u_i

A=@(i) dt.*[1/dt 0 -U_ref(1,i)*sin(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*cos(Y_ref(3,i));
            0 1/dt  U_ref(1,i)*cos(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*sin(Y_ref(3,i));
            0 0 1/dt];

B=@(i) dt.*[cos(Y_ref(3,i))-b/L*tan(U_ref(2,i))*sin(Y_ref(3,i)) -b/L*U_ref(1,i)*sin(Y_ref(3,i))/cos(U_ref(2,i))^2;
            sin(Y_ref(3,i))+b/L*tan(U_ref(2,i))*cos(Y_ref(3,i)) b/L*U_ref(1,i)*cos(Y_ref(3,i))/cos(U_ref(2,i))^2;
            tan(U_ref(2,i))/L U_ref(1,i)/L/cos(U_ref(2,i))^2];


%% 2.2 Number of decision variables for colocation method
npred=10;

Ndec= 11*3 + 10*2;

%% 2.3 write and test function  to construct Aeq beq (equality constraints
%enforce x(i+1)=A_i*x_i+B_i*u_i for the prediction horizon) 
%check the values of Aeq and beq timestep 1 (t=0)
%with the given initial condition of Y(0)=[0.25;-0.25;-0.1];

[Aeq,beq]=eq_cons(1,A,B,10, [10,10,pi]);
Aeq_test1= Aeq;
beq_test1= beq;

%% 2.4 write and test function to generate limits on inputs 
%check the value at timestep 1 (t=0)
[Lb_test1,Ub_test1]=bound_cons(1,U_ref,10);

%% 2.5 simulate controller working from initial condition [0.25;-0.25;-0.1]
%use ode45 to between inputs
Q = [1 0 0; 0 1 0; 0 0 1];
R = [0.1 0; 0 0.1];

x = zeros(3, length(T));
u = zeros(length(T)-1 , 2);
x(:, 1) = X0;
PredHorizon = 11;

H = zeros(3* PredHorizon + 2*(PredHorizon - 1) );
c = zeros(3* PredHorizon + 2*(PredHorizon - 1), 1);
for i = 1:PredHorizon
    H(((i-1)*3+1):(i*3), ((i-1)*3+1):(i*3)) = Q;
end
for i = 1:(PredHorizon - 1)
    H( ((i-1)*2 + 3*PredHorizon + 1 ):( i*2 + 3*PredHorizon), ((i-1)*2 + 3*PredHorizon+1 ):(i*2 + 3*PredHorizon) ) = R;
end


odefun = @(t, x, u) [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
                  u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
                  u(1)*tan(u(2))/L] ;
 
  odefun_d=@(t,x, u) [
    x(4)*cos(x(3))-x(5)*sin(x(3));
    x(4)*sin(x(3))+x(5)*cos(x(3));
    x(6);
    1/m*((Cm1*u(2)-Cm2*u(2)*x(4)-Cr0-Cr2*x(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*sin(u(1)) + m*x(5)*x(6));
 ( 1/m*(D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))) + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*cos(u(1)) - m*x(4)*x(6))) ;
( 1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*l_f*cos(u(1))- (D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))))*l_r)  )
    ];
              
              
j = 1;

%% parameter setup for Kalman
X0_d_m=X0_d;

N_R=[0.49 0.49 0.01 0.01 0.0001 0.0001];
N_Q = diag([0 0 0 0.01*0.01 0.01*0.01 0.01*0.01]);
noise=zeros(6, length(T));
xhat= [X0_d,X0_d,X0_d,X0_d,X0_d,X0_d];
xhat_stor=zeros(6,length(T));
xhat_stor(:,1)=X0_d;

%%


while(j < (length(T))-100)
    if(npred > (length(T)-j))
        npred = min(npred, (length(T)-j));
        PredHorizon = npred+1;
        H = zeros(3* PredHorizon + 2*(PredHorizon - 1) );
        c = zeros(3* PredHorizon + 2*(PredHorizon - 1), 1);
        for i = 1:PredHorizon
            H(((i-1)*3+1):(i*3), ((i-1)*3+1):(i*3)) = Q;
        end
        for i = 1:(PredHorizon - 1)
            H( ((i-1)*2 + 3*PredHorizon + 1 ):( i*2 + 3*PredHorizon), ((i-1)*2 + 3*PredHorizon+1 ):(i*2 + 3*PredHorizon) ) = R;
        end
    end 
    
%%
%Kalman fitler in loop
   
if j>1
        
for p=1:1:6
noise(p,j)= X0_d(p)+normrnd(0,sqrt(N_R(p)));
end
    
    if Kalman_on==1

    sys_dis=linearization(xhat_stor(:,j-1), u_dyna, dt);

    output(:,1)=ExtKalman1(noise(1,j), sys_dis,N_R,N_Q,1,xhat(:,1),u_dyna,dt);
    output(:,2)=ExtKalman2(noise(2,j), sys_dis,N_R,N_Q,2,xhat(:,2),u_dyna,dt);
    output(:,3)=ExtKalman3(noise(3,j), sys_dis,N_R,N_Q,3,xhat(:,3),u_dyna,dt);
    output(:,4)=ExtKalman4(noise(4,j), sys_dis,N_R,N_Q,4,xhat(:,4),u_dyna,dt);
    output(:,5)=ExtKalman5(noise(5,j), sys_dis,N_R,N_Q,5,xhat(:,5),u_dyna,dt);
    output(:,6)=ExtKalman6(noise(6,j), sys_dis,N_R,N_Q,6,xhat(:,6),u_dyna,dt);
    xhat=[output(:,1),output(:,2),output(:,3),output(:,4),output(:,5),output(:,6)] ;
    X0_d_m=[xhat(1,1);xhat(2,2) ;xhat(3,3) ;xhat(4,4) ;xhat(5,5); xhat(6,6)];
    xhat_stor(:,j)=X0_d_m;



%plot(stateseries(1:step+1,1),stateseries(1:step+1,2),'*');

%axis([-5,10,0,10])

%pause(0.1)
    else
      X0_d_m= noise(:,j);
    end    

        
end
    %%
    [Aeq, beq] = eq_cons(j, A, B, npred, X0_d_m(1:3)-Y_ref(:, j)); 
    [Lb,Ub]=bound_cons(j, U_ref, npred);
    solDV = quadprog( H, c, [], [], Aeq, beq, Lb, Ub);
    u = solDV((3*PredHorizon + 1):(2 + 3* PredHorizon))+ U_ref(:, j);
    u(1)-X0_d(4);
    u_d=Inversesolver_d(u,X0_d, dt);
    u_dyna=[u(2) u_d];

    [~, temp_x] = ode45(@(t,x) odefun(t, x, u), [0 dt], X0);
    [~, temp_xd] = ode45(@(t,x) odefun_d(t, x, u_dyna), [0 dt], X0_d);
    X0 = temp_x(end, :)'; 
    X0_d= temp_xd(end, :)'; 
    x(:, j+1) = solDV(4:6) + Y_ref(:, j+1);
    x_store(:,j+1)=X0;
    x_store_d(:,j+1)=X0_d;
    j = j + 1


end

max_dist_error = 0;
for t = 2:1:length(T)  
    if(x(1, t) >= 3 && x(1, t) <= 4)
        cur_error = sqrt((x(1, t)- Y_ref(1, t)).^2 + (x(2, t)- Y_ref(2, t)).^2);
        max_dist_error = max(max_dist_error, cur_error);
    end
end
%plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); hold on;
%plot(x(1, 1:length(x) - 10), x(2, 1:length(x) - 10),'r'); hold on;
%plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); hold on;
%plot(x_store(1, 1:length(x_store) - 10), x_store(2, 1:length(x_store) - 10),'r'); hold on;

if Kalman_on==1
    hold on;
    plot(x_store_d(1, 1:1600), x_store_d(2, 1:1600),'b--');
plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'k');
legend('Noise with fitler','reference')
else
    hold on;
    plot(x_store_d(1, 1:1600), x_store_d(2, 1:1600),'r--'); 
plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'k');
legend('Noise without fitler','reference')
end
    
hold on;




%Obstacle1_x = [1,3,3,1,1,3,1,3];
%Obstacle1_y = [-1,-1,1,1,-1,1,1,-1]*2.25;
%plot(Obstacle1_x, Obstacle1_y,'k');
%hold off;
%

%% I highly recommend you write functions for constraint generation and dynamics down here and 
%call them when needed above, for exa

function [Aeq,beq]=eq_cons(idx,A,B,step, IC)
Aeq = zeros(33,53);
Aeq(1,1) = 1;
Aeq(2,2) = 1;
Aeq(3,3) = 1;
beq = zeros(33,1);
beq(1) = IC(1);
beq(2) = IC(2);
beq(3) = IC(3);

for t = 2:step+1
    A_temp = A(idx + t - 2);
    B_temp = B(idx + t - 2);
    Aeq(3*t-2,3*t-2) = -1;
    Aeq(3*t-2,3*t-5) = A_temp(1,1);
    Aeq(3*t-2,3*t-4) = A_temp(1,2);
    Aeq(3*t-2,3*t-3) = A_temp(1,3);
    Aeq(3*t-2,30+2*t) = B_temp(1,1);
    Aeq(3*t-2,31+2*t) = B_temp(1,2);
    
    Aeq(3*t-1,3*t-1) = -1;   
    Aeq(3*t-1,3*t-5) = A_temp(2,1);
    Aeq(3*t-1,3*t-4) = A_temp(2,2);
    Aeq(3*t-1,3*t-3) = A_temp(2,3);
    Aeq(3*t-1,30+2*t) = B_temp(2,1);
    Aeq(3*t-1,31+2*t) = B_temp(2,2);
    
    Aeq(3*t,3*t) = -1;   
    Aeq(3*t,3*t-5) = A_temp(3,1);
    Aeq(3*t,3*t-4) = A_temp(3,2);
    Aeq(3*t,3*t-3) = A_temp(3,3);
    Aeq(3*t,30+2*t) = B_temp(3,1);
    Aeq(3*t,31+2*t) = B_temp(3,2);
    
end
%build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
%in the form Aeq*z=beq
%initial_idx specifies the time index of initial condition from the reference trajectory 
%A and B are function handles above
end

%function [Lb,Ub]=bound_cons(idx,U_ref,input_range)
%initial_idx is the index along uref the initial condition is at
function [Lb,Ub]=bound_cons(idx,U_ref,step)
    Lb = zeros(53,1);
    Ub = zeros(53,1);
    for t = 1:1+step
        Lb(3*t-2) = -Inf;
        Ub(3*t-2) = Inf;
        Lb(3*t-1) = -Inf;
        Ub(3*t-1) = Inf;
        Lb(3*t) = -Inf;
        Ub(3*t) = Inf;
    end
    for t = 1:step
        Lb(32+2*t) = 0-U_ref(1,idx+t-1);
        Ub(32+2*t) = 15-U_ref(1,idx+t-1);
        Lb(33+2*t) = -1-U_ref(2,idx+t-1);
        Ub(33+2*t) = 1-U_ref(2,idx+t-1);
    end

end
