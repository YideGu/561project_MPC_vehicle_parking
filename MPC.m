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
T=0:dt:18;

X0 = [7.1 7.5 pi]';

%% load reference trajectory
k = 1.15;
h = k/300;
U1 = [0:h:k];
U2 = [k-h:-h:0];
%U3 = [k-h:-h:h, 0.002];
Ut = [U1, U2];

%{
U_ref(1:2,1:601)=[ones(1,601) * 1; zeros(1,601)];

[x,y,phi] = kinematic_model(6,[9.2 10 pi]',b,L,U_ref(1:2,1:601));

Y_ref(1:3,1:601)=[x;y;phi];

U_ref(1:2,601:1201)=[ones(1,601) * 1.5; Ut];

[x,y,phi] = kinematic_model(6,Y_ref(:,601),b,L,U_ref(1:2,601:1201));

Y_ref(1:3,601:1201)=[x;y;phi];

U_ref(1:2,1201:1801)=[ones(1,601) * 0.3; ones(1,601) * 0]; %0.6

[x,y,phi] = kinematic_model(6,Y_ref(:,1201),b,L,U_ref(1:2,1201:1801));

Y_ref(1:3,1201:1801)=[x;y;phi];
%}
% U_ref(1:2,1:601)=[ones(1,601) * 1; zeros(1,601)];
% 
% [x,y,phi] = kinematic_model(6,[6.5 10.2 pi]',b,L,U_ref(1:2,1:601));
% 
% Y_ref(1:3,1:601)=[x;y;phi];
% 
% U_ref(1:2,601:1201)=[ones(1,601) * 0.675; ones(1,601) * 1.05];
% 
% [x,y,phi] = kinematic_model(6,Y_ref(:,601),b,L,U_ref(1:2,601:1201));
% 
% Y_ref(1:3,601:1201)=[x;y;phi];
% 
% U_ref(1:2,1201:1801)=[ones(1,601) * 1; ones(1,601) * 0]; %0.6
% 
% [x,y,phi] = kinematic_model(6,Y_ref(:,1201),b,L,U_ref(1:2,1201:1801));
% 
% Y_ref(1:3,1201:1801)=[x;y;phi];
% 
% plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); 
% hold on;
U_ref = zeros(2,1801);
U_ref(1:2,1:401) = [ones(1,401)*1;ones(1,401)*0];
U_ref(1:2,401:1201) = [ones(1,801)* 1.048; ones(1,801) *0.7];
U_ref(1:2,1201:1801) = [ones(1,601)* 0.8; ones(1,601) * 0];

[x,y,phi] = kinematic_model(18,[7.1 12.38 pi]',b,L,U_ref,dt);
Y_ref=[x;y;phi];
plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); 


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
Q = [1 0 0; 0 1 0; 0 0 0.5];
R = [0.1 0; 0 0.01];

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
 
j = 1;
while(j < (length(T))-10)
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
    
    [Aeq, beq] = eq_cons(j, A, B, npred, X0-Y_ref(:, j)); 
    [Lb,Ub]=bound_cons(j, U_ref, npred);
    solDV = quadprog( H, c, [], [], Aeq, beq, Lb, Ub);
    u = solDV((3*PredHorizon + 1):(2 + 3* PredHorizon))+ U_ref(:, j);
    [~, temp_x] = ode45(@(t,x) odefun(t, x, u), [0 dt], X0);
    X0 = temp_x(end, :)'; 
    x(:, j+1) = solDV(4:6) + Y_ref(:, j+1);
    j = j + 1; 
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
        Lb(33+2*t) = -1.0-U_ref(2,idx+t-1);
        Ub(33+2*t) = 1.0-U_ref(2,idx+t-1);
    end

end
