function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Your code here
Y = zeros(length(R),size(X,2)-1);                   
a = (mvnrnd(zeros(length(R),1),R,size(X,2)-1))';       
for i = 1:(size(X,2)-1)
    [hx,~] = h(X(:,i+1));                       
    Y(:,i) = hx+a(:,i);
end
end