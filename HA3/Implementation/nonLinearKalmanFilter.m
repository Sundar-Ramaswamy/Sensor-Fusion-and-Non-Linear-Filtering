function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%                       [fx,Fx]=f(x) 
%                       Takes as input x (state) 
%                       Returns fx and Fx, motion model and Jacobian evaluated at x
%   Q           [n x n] Process noise covariance
%   h                   Measurement model function handle
%                       [hx,Hx]=h(x,T) 
%                       Takes as input x (state), 
%                       Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

% Your code here. If you have good code for the Kalman filter, you should re-use it here as
% much as possible.
N = size(Y,2);
%% Data allocation
xp = zeros(length(x_0),N);
Pp = zeros(length(x_0),length(x_0),N);
xf = zeros(length(x_0),N+1);
Pf = zeros(length(x_0),length(x_0),N+1);
%% Filter Implementation
xf(:,1)=x_0;
Pf(:,:,1)=P_0;
for i = 1:N
    [xp(:,i), Pp(:,:,i)] = nonLinKFprediction(xf(:,i), Pf(:,:,i), f, Q, type);
    [xf(:,i+1), Pf(:,:,i+1)] = nonLinKFupdate(xp(:,i), Pp(:,:,i), Y(:,i), h, R, type);
end    
xf = xf(:,2:end);
Pf = Pf(:,:,2:end);
end