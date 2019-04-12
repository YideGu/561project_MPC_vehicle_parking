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
T=0.05;                         % discrete time
state=[10,10,pi, 3, 0, 0];      % initial state
num=40;                         % the time of the calculation
stateseries=zeros(num+1,6);     % the vector to store the state over time
stateseries(1,:)=state;

% plot designed trajectory 
hold on 
xline = 0: pi/100:3*pi;
yline = 10*sin(xline/6);
plot(xline,yline,'--')

%plot actual trajectory 
for step=1:1:num
%%%%%%%%% this input needs designed
u=[0.2+0.03*step,0.6];
%%%%%%%%%%%%

[t,z]=trajectory(state,u,pf,p,T);   % return the next time state associated with exact nonlinear system 
stateseries(step+1,:)=z(end,:);

plot(stateseries(1:step+1,1),stateseries(1:step+1,2),'*');
axis([-5,10,0,10])
state=z(end,:);
pause(0.1)
end
