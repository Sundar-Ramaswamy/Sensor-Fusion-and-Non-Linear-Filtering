function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% your code here
Y = zeros(size(H,1),size(X,2)-1);                       
r = (mvnrnd(zeros(size(R,1),1),R,size(X,2)-1))';       
for i = 1:(size(X,2)-1)
    Y(:,i) = H*X(:,i+1)+r(:,i);   
end
end