function xhatOut = ExtKalman5(meas,sys_d,R,Q,i,xhat,u,T)
% This Embedded MATLAB Function implements an extended Kalman filter used
% for object tracking.
%
% The states of the process are given by
% x = [x_position; x_velocity; y_position; y_velocity];
%
% and the measurements are given by
% y = [range; bearing]
%
% where
% range = sqrt(x_position^2 + y_position^2)
% bearing = atan2(y_position/x_position)

% Author: Phil Goddard (phil@goddardconsulting.ca)
% Date: Q2, 2011.

% Define storage for the variables that need to persist
% between time periods.


persistent P5 
   
  
if isempty(P5)  % First time through the code so do some initialization
 
   P5 = zeros(6,6);
  
   
end
   A=sys_d.A;
   B=sys_d.B;
% Calculate the Jacobians for the state and measurement equations
yhat=xhat(i);
%
        H=[0,0,0,0,1,0];
% Propogate the state and covariance matrices
xhat = xhat+T*differe(xhat,u(1),u(2));
P5 = A*P5*A' + Q;
% Calculate the Kalman gain
K = P5*H'/(H*P5*H' + R(i));
% Calculate the measurement residual
resid = meas - yhat;
% Update the state and covariance estimates
xhat = xhat + K*resid;
P5 = (eye(size(K,1))-K*H)*P5;
% Post the results

xhatOut = xhat;