function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
%jointGaussian calculates the joint Gaussian density as defined
%in problem 1.3a. 
%
%Input
%   MU_X        Expected value of x
%   SIGMA2_X    Covariance of x
%   SIGMA2_R    Covariance of the noise r
%
%Output
%   MU          Mean of joint density 
%   SIGMA       Covariance of joint density


%Your code here
A1=[1 0; 1 1];
b1=[0;0];
mu=A1*[mu_x; 0]+b1
Sigma=A1*blkdiag(sigma2_x,sigma2_r)*A1'


end