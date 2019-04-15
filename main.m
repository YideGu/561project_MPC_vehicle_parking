%%%%%%%%%%%%%%%%%%%%%%%
%WI19 VM 561 semester project %
%This file is the mian file to run 
%
%
%
%%
clear all

% initialize the vehicle-related parameter from another file
parametersetup;

% determine algorithm parameter
T=0.002;                         % discrete time
state=[12,10,pi, 3, 0, 0.5];      % initial state
num=500;                         % the time of the calculation
stateseries=zeros(num+1,6);     % the vector to store the state over time
stateseries(1,:)=state;

% plot designed trajectory 
hold on 
xline = 0: pi/100:3*pi;
yline = 10*sin(xline/6);
plot(xline,yline,'--')
raw=zeros(6,num);
meas_stor=zeros(6,num);
xhat_stor=zeros(6,num);
output=zeros(6,6);
%plot actual trajectory 
R=[1 1 0.01 0.01 0.001 0.001];
Q = diag([0 0 0 0.001*0.001 0.001*0.001 0.001*0.001]);

xhat= [state',state',state',state',state',state'];
for step=1:1:num
%%%%%%%%% this input needs designed
u=[0.2+0.3/num*step,0.6-0.5/num*step];
%%%%%%%%%%%%
sys_dis=linearization(state, u, pf,p, T);


[t,z]=trajectory(state,u,pf,p,T);   % return the next time state associated with exact nonlinear system 

for j=1:1:6
raw(j,step)= z(end,j);
meas_stor(j,step)= raw(j,step)+normrnd(0,sqrt(R(j)));
end
output(:,1)=ExtKalman1(meas_stor(1,step), sys_dis,R,Q,1,xhat(:,1),u',T);
output(:,2)=ExtKalman2(meas_stor(2,step), sys_dis,R,Q,2,xhat(:,2),u',T);
output(:,3)=ExtKalman3(meas_stor(3,step), sys_dis,R,Q,3,xhat(:,3),u',T);
output(:,4)=ExtKalman4(meas_stor(4,step), sys_dis,R,Q,4,xhat(:,4),u',T);
output(:,5)=ExtKalman5(meas_stor(5,step), sys_dis,R,Q,5,xhat(:,5),u',T);
output(:,6)=ExtKalman6(meas_stor(6,step), sys_dis,R,Q,6,xhat(:,6),u',T);
xhat=[output(:,1),output(:,2),output(:,3),output(:,4),output(:,5),output(:,6)] ;
xhat_stor(:,step)=[xhat(1,1);xhat(2,2) ;xhat(3,3) ;xhat(4,4) ;xhat(5,5); xhat(6,6)];

stateseries(step+1,:)=z(end,:);


%plot(stateseries(1:step+1,1),stateseries(1:step+1,2),'*');

%axis([-5,10,0,10])
state=z(end,:);
%pause(0.1)
end
xaxis= [1:1:step];

for k=1:1:6
    subplot(6,1,k)
plot(xaxis,raw(k,:),xaxis,xhat_stor(k,:),'--',xaxis,meas_stor(k,:))
xlim([0,500])
legend('real','estimate','measure')
end
