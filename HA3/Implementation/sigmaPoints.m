function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%
n = length(x);
    switch type        
        case 'UKF'
    
            SP = zeros(n,2*n+1);                                                                                      
            SP(:,1) = x;                                          
            Sqrt_P = sqrtm(P);
            W0 = 1-n/3;
            for i = 1:size(P,2)
                SP(:,i+1) = x + sqrt((n)/(1-W0))*Sqrt_P(:,i);        
                SP(:,i+n+1) = x - sqrt((n)/(1-W0))*Sqrt_P(:,i);     
            end
            W = [W0, ((1-W0)/(2*n))*ones(1,2*n)];
                
        case 'CKF'
            
            SP = zeros(n,2*n);                                      
            Sqrt_P = sqrtm(P);
            for i = 1:size(P,2)
                SP(:,i) = x + sqrt(n)*Sqrt_P(:,i);                   
                SP(:,i+n) = x - sqrt(n)*Sqrt_P(:,i);                
            end
            W = ((1)/(2*n))*ones(1,2*n);
            
        otherwise
            error('Incorrect type of sigma point')
    end

end