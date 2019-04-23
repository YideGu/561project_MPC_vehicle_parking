%%This script is for Experiment 1
%Apply MPC on kinematic model, and get trajectory using kinematic model
clear;
clc;
%% parameters
%wheelbases
%distance between front and rear wheel 
L=2.7;
%distance from rear wheel to center of mass
b=1.35;

%time discretization
dt=0.01;
%time span We make the simulation in 18 seconds
T=0:dt:18;

X0 = [4.1 7.5 pi]';

%% load reference trajectory
% Set a pre-given input reference 
U_ref = zeros(2,1801);
U_ref(1:2,1:401) = [ones(1,401)*1;ones(1,401)*0];
U_ref(1:2,401:1201) = [ones(1,801)* 0.63; ones(1,801) *0.7];
U_ref(1:2,1201:1801) = [ones(1,601)* 0.8; ones(1,601) * 0];

%Use kinematic model to compute state reference from input reference
[x,y,phi] = kinematic_model(18,[5.84 9.36 pi]',b,L,U_ref);
Y_ref=[x;y;phi];
plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); 


%% Discrete-time A and B matrices
%these are the system linearized in discrete time about the reference
A=@(i) dt.*[1/dt 0 -U_ref(1,i)*sin(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*cos(Y_ref(3,i));
            0 1/dt  U_ref(1,i)*cos(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*sin(Y_ref(3,i));
            0 0 1/dt];

B=@(i) dt.*[cos(Y_ref(3,i))-b/L*tan(U_ref(2,i))*sin(Y_ref(3,i)) -b/L*U_ref(1,i)*sin(Y_ref(3,i))/cos(U_ref(2,i))^2;
            sin(Y_ref(3,i))+b/L*tan(U_ref(2,i))*cos(Y_ref(3,i)) b/L*U_ref(1,i)*cos(Y_ref(3,i))/cos(U_ref(2,i))^2;
            tan(U_ref(2,i))/L U_ref(1,i)/L/cos(U_ref(2,i))^2];


%% Number of decision variables for colocation method
npred=10;

Ndec= 11*3 + 10*2;

%% simulate controller working from initial condition
%use ode45 to between inputs
Q = [1 0 0; 0 1 0; 0 0 0.5];
R = [0.1 0; 0 0.01];

x = zeros(3, length(T));
u = zeros(length(T)-1 , 2);
x(:, 1) = X0;
PredHorizon = 11;

odefun = @(t, x, u) [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
                  u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
                  u(1)*tan(u(2))/L] ;
              
H = zeros(3* PredHorizon + 2*(PredHorizon - 1) );
c = zeros(3* PredHorizon + 2*(PredHorizon - 1), 1);
for i = 1:PredHorizon
    H(((i-1)*3+1):(i*3), ((i-1)*3+1):(i*3)) = Q;
end
for i = 1:(PredHorizon - 1)
    H( ((i-1)*2 + 3*PredHorizon + 1 ):( i*2 + 3*PredHorizon), ((i-1)*2 + 3*PredHorizon+1 ):(i*2 + 3*PredHorizon) ) = R;
end
 
              
              
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


Obstacle1_x = [1.3,3.9,3.9,1.3,1.3,3.9,1.3,3.9];
Obstacle1_y = [-1,-1,1,1,-1,1,1,-1] * 2.7;
Obstacle2_x = [1.3,3.9,3.9,1.3,1.3,3.9,1.3,3.9] * -1;
plot(Obstacle1_x, Obstacle1_y,'k');hold on; 
plot(Obstacle2_x, Obstacle1_y,'k');hold on;

plot(Y_ref(1, 1: length(Y_ref) - 10), Y_ref(2, 1:length(Y_ref) - 10),'b'); hold on;
plot(x(1, 1:length(x) - 10), x(2, 1:length(x) - 10),'r'); hold on;
legend('reference','kinematic trajectory')

hold off;


%% Error dynamics
max_dist_error = 0;
for t = 2:1:length(T)  
    if(x(2, t) <= 6)
        cur_error = sqrt((x(1, t)- Y_ref(1, t)).^2 + (x(2, t)- Y_ref(2, t)).^2);
        max_dist_error = max(max_dist_error, cur_error);
    end
end
disp(max_dist_error)

%% Self defined function
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


end


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
