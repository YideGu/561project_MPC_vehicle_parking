clear;
clc;
%% parameters
%wheelbases
% bounding box 4.5m x 1.7m for compact car

L=4.5;
%distance from rear wheel to center of mass
b=2.25;

%time discretization
dt=0.01;
%time span We make the simulation in 18 seconds
T=0:dt:6;

X0 = [6.5 10.2 0 3 0 0]';

%% load reference trajectory


K=0.01;
K_n=100;
G=0.6;
G_n=187;

Y_ref(:,1)=[6.5 10.2 0 3 0 0];

U_ref(1:2,1:200)=[zeros(1,200); [K*ones(1,K_n), zeros(1, 200-K_n)]];

[x,y,phi, v_x,v_y,omega] = dynamic_model(dt,Y_ref(:,1)',U_ref(1:2,1:200));

Y_ref(1:6,2:201)=[x;y;phi; v_x;v_y;omega];

U_ref(1:2,201:401)=[[G*ones(1,G_n), zeros(1,201-G_n)];  0.1*K*ones(1,201) ];

[x,y,phi, v_x,v_y,omega] = dynamic_model(dt,Y_ref(:,201),U_ref(1:2,201:401));

Y_ref(1:6,201:401)=[x;y;phi; v_x;v_y;omega];

U_ref(1:2,401:601)=[ zeros(1,201);[K*ones(1,K_n), zeros(1, 201-K_n)]]; %0.6

[x,y,phi, v_x,v_y,omega] = dynamic_model(dt,Y_ref(:,401),U_ref(1:2,401:601));

Y_ref(1:6,401:601)=[x;y;phi; v_x;v_y;omega];
hold on

plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); 
plot(Y_ref(1, 1: 400), Y_ref(2, 1:400),'r'); 
plot(Y_ref(1, 1: 200), Y_ref(2, 1:200),'y'); 




Animation(Y_ref)
%%



%the reference trajectory Yref is given in a 3x601 double

%% 2.1 Discrete-time A and B matrices
%these are the system linearized in discrete time about the reference
%trajectory i.e. x(i+1)=A_i*x_i+B_i*u_i
parametersetup
[xxx,Afun,Bfun]=linearization([10,10,0,1, 0, 0],[1,0],pf,p,dt)

%{
A=@(i) dt.*[1/dt 0 -U_ref(1,i)*sin(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*cos(Y_ref(3,i));
            0 1/dt  U_ref(1,i)*cos(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*sin(Y_ref(3,i));
            0 0 1/dt];

B=@(i) dt.*[cos(Y_ref(3,i))-b/L*tan(U_ref(2,i))*sin(Y_ref(3,i)) -b/L*U_ref(1,i)*sin(Y_ref(3,i))/cos(U_ref(2,i))^2;
            sin(Y_ref(3,i))+b/L*tan(U_ref(2,i))*cos(Y_ref(3,i)) b/L*U_ref(1,i)*cos(Y_ref(3,i))/cos(U_ref(2,i))^2;
            tan(U_ref(2,i))/L U_ref(1,i)/L/cos(U_ref(2,i))^2];
%}
linear_dynamic
        
        
%% 2.2 Number of decision variables for colocation method
npred=10;

Ndec= 11*3 + 10*2;

%% 2.3 write and test function  to construct Aeq beq (equality constraints
%enforce x(i+1)=A_i*x_i+B_i*u_i for the prediction horizon) 
%check the values of Aeq and beq timestep 1 (t=0)
%with the given initial condition of Y(0)=[0.25;-0.25;-0.1];

[Aeq,beq]=eq_cons(1,A,B,10, [10,10,0,0, 0, 0]);
Aeq_test1= Aeq;
beq_test1= beq;

%% 2.4 write and test function to generate limits on inputs 
%check the value at timestep 1 (t=0)
[Lb_test1,Ub_test1]=bound_cons(1,U_ref,10);

%% 2.5 simulate controller working from initial condition [0.25;-0.25;-0.1]
%use ode45 to between inputs
Q = [1 0 0 0 0 0; 0 1 0 0 0 0 ; 0 0 0.1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0;0 0 0 0 0 0.1];
R = [0.5 0; 0 0.5];

x = zeros(6, length(T));
u = zeros(length(T)-1 , 2);
x(:, 1) = X0;
PredHorizon = 11;

H = zeros(6* PredHorizon + 2*(PredHorizon - 1) );
c = zeros(6* PredHorizon + 2*(PredHorizon - 1), 1);
for i = 1:PredHorizon
    H(((i-1)*6+1):(i*6), ((i-1)*6+1):(i*6)) = Q;
end
for i = 1:(PredHorizon - 1)
    H( ((i-1)*2 + 6*PredHorizon + 1 ):( i*2 + 6*PredHorizon), ((i-1)*2 + 6*PredHorizon+1 ):(i*2 + 6*PredHorizon) ) = R;
end


%odefun = @(t, x, u) [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
%                  u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
%                  u(1)*tan(u(2))/L] ;
 
       
 
odefun=@(t,x, u) [
    x(4)*cos(x(3))-x(5)*sin(x(3));
    x(4)*sin(x(3))+x(5)*cos(x(3));
    x(6);
    1/m*((Cm1*u(2)-Cm2*u(2)*x(4)-Cr0-Cr2*x(4)^2) - (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*sin(u(1)) + m*x(5)*x(6));
sign(u(1))* ( 1/m*(D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))) + (D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*cos(u(1)) - m*x(4)*x(6)) );
sign(u(1))* ( 1/Iz*((D_f*sin(C_f*atan(B_f*( -atan2(l_f*x(6) + x(5),x(4))+u(1)))))*l_f*cos(u(1))- (D_r*sin(C_r*atan(B_r*(atan2(l_r*x(6) - x(5),x(4))))))*l_r)  )
    ];

              
j = 1;

while(j < (length(T))-10)

    if(npred > (length(T)-j))
        npred = min(npred, (length(T)-j));
        PredHorizon = npred+1;
        H = zeros(6* PredHorizon + 2*(PredHorizon - 1) );
        c = zeros(6* PredHorizon + 2*(PredHorizon - 1), 1);
        for i = 1:PredHorizon
            H(((i-1)*6+1):(i*6), ((i-1)*6+1):(i*6)) = Q;
        end
        for i = 1:(PredHorizon - 1)
            H( ((i-1)*2 + 6*PredHorizon + 1 ):( i*2 + 6*PredHorizon), ((i-1)*2 + 6*PredHorizon+1 ):(i*2 + 6*PredHorizon) ) = R;
        end
    end
    
    [Aeq, beq] = eq_cons(j, A, B, npred, X0-Y_ref(:, j)); 
    devia(:,j)=X0-Y_ref(:, j);
    [Lb,Ub]=bound_cons(j, U_ref, npred);
    solDV = quadprog( H, c, [], [], Aeq, beq, Lb, Ub);
    u = solDV((6*PredHorizon + 1):(2 + 6* PredHorizon))+ U_ref(:, j)
    ustore(:,j)=u;
    [~, temp_x] = ode23(@(t,x) odefun(t, x, u), [0 dt], X0);
    X0 = temp_x(end, :)';
    
    x(:, j+1) = solDV(7:12) + Y_ref(:, j+1);
    j = j + 1
end


max_dist_error = 0;
for t = 2:1:length(T)  
    if(x(1, t) >= 3 && x(1, t) <= 4)
        cur_error = sqrt((x(1, t)- Y_ref(1, t)).^2 + (x(2, t)- Y_ref(2, t)).^2);
        max_dist_error = max(max_dist_error, cur_error);
    end
end
plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); hold on;
plot(x(1, 1:length(x) - 10), x(2, 1:length(x) - 10),'r'); hold on;

Obstacle1_x = [1,3,3,1,1,3,1,3];
Obstacle1_y = [-1,-1,1,1,-1,1,1,-1]*2.25;
plot(Obstacle1_x, Obstacle1_y,'k');
hold off;

%% I highly recommend you write functions for constraint generation and dynamics down here and 
%call them when needed above, for exa

function [Aeq,beq]=eq_cons(idx,A,B,step, IC)
Aeq = zeros(66,86);
Aeq(1,1) = 1;
Aeq(2,2) = 1;
Aeq(3,3) = 1;
Aeq(4,4) = 1;
Aeq(5,5) = 1;
Aeq(6,6) = 1;
beq = zeros(66,1);
beq(1) = IC(1);
beq(2) = IC(2);
beq(3) = IC(3);
beq(4) = IC(4);
beq(5) = IC(5);
beq(6) = IC(6);
for t = 2:step+1
    A_temp = A(idx + t - 2);
    B_temp = B(idx + t - 2);
    for j=0:1:5
    Aeq(6*t-j,6*t-j) = -1;
    Aeq(6*t-j,6*t-11) = A_temp(6-j,1);
    Aeq(6*t-j,6*t-10) = A_temp(6-j,2);
    Aeq(6*t-j,6*t-9) = A_temp(6-j,3);
    Aeq(6*t-j,6*t-8) = A_temp(6-j,4);
    Aeq(6*t-j,6*t-7) = A_temp(6-j,5);
    Aeq(6*t-j,6*t-6) = A_temp(6-j,6);
    Aeq(6*t-j,63+2*t) = B_temp(6-j,1);
    Aeq(6*t-j,64+2*t) = B_temp(6-j,2);
    end
%{    
    Aeq(3*t,3*t) = -1;   
    Aeq(3*t,3*t-5) = A_temp(3,1);
    Aeq(3*t,3*t-4) = A_temp(3,2);
    Aeq(3*t,3*t-3) = A_temp(3,3);
    Aeq(3*t,30+2*t) = B_temp(3,1);
    Aeq(3*t,31+2*t) = B_temp(3,2);
%}    
end
%build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
%in the form Aeq*z=beq
%initial_idx specifies the time index of initial condition from the reference trajectory 
%A and B are function handles above


end

%function [Lb,Ub]=bound_cons(idx,U_ref,input_range)
%initial_idx is the index along uref the initial condition is at
function [Lb,Ub]=bound_cons(idx,U_ref,step)
    Lb = zeros(86,1);
    Ub = zeros(86,1);
    for t = 1:1+step
        Lb(6*t-5) = -Inf;
        Ub(6*t-5) = Inf;
        Lb(6*t-4) = -Inf;
        Ub(6*t-4) = Inf;
        Lb(6*t-3) = -Inf;
        Ub(6*t-3) = Inf;
        Lb(6*t-2) = -Inf;
        Ub(6*t-2) = Inf;
        Lb(6*t-1) = -Inf;
        Ub(6*t-1) = Inf;
        Lb(6*t) = -Inf;
        Ub(6*t) = Inf;
    end
    
    
    for t = 1:step
        Lb(65+2*t) = -1-U_ref(1,idx+t-1);
        Ub(65+2*t) = 1-U_ref(1,idx+t-1);
        Lb(66+2*t) = -1-U_ref(2,idx+t-1);
       Ub(66+2*t) = 1-U_ref(2,idx+t-1);
   end
%{
    for t = 1:step
        Lb(65+2*t) = -Inf;
        Ub(65+2*t) = Inf;
        Lb(66+2*t) = -Inf;
       Ub(66+2*t) = Inf;
    end
%}
    
end