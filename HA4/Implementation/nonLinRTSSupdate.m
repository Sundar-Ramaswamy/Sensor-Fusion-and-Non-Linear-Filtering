function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

% Your code here.
switch type
 case 'EKF'
 [fx, Fx]=f(xf_k,T);
 P_k=Pf_k*Fx';
 update = xs_kplus1-xp_kplus1; 
 case 'UKF'
 [SP,W]=sigmaPoints(xf_k,Pf_k,type); P_k=0;
 for i=1:length(SP)
 P_k=P_k+(SP(:,i)-xf_k)*(f(SP(:,i),T)-xp_kplus1)'*W(i);
 end
 update=xs_kplus1-xp_kplus1;
 case 'CKF'
 [SP,W]=sigmaPoints(xf_k, Pf_k,type);
 P_k=0;
 
 for i=1:length(SP)
     P_k=P_k+(SP(:,i)-xf_k)*(f(SP(:,i),T)-xp_kplus1).'*W(i);
 end
 update=xs_kplus1-xp_kplus1;
end
G_k = P_k * inv(Pp_kplus1);
xs = xf_k + G_k * update;
Ps = Pf_k - G_k * (Pp_kplus1 - Ps_kplus1) * G_k';
end
