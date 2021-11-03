function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Non-resampled Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
n=size(x_0,1); [m,K]=size(Y);
Xp=zeros(n,N,K); Wp=zeros(N,K); 
if nargout > 2 
    xfp=zeros(n,K); Pfp=zeros(n,n,K);
end
%Non-resampled weights for the posterior state x
Wp(:,1)=repmat(1/N,[N,1]);
%Particle for the posterior state from prior
Xp(:,:,1)=mvnrnd(x_0,P_0,N)';
%Particle for the posterior state from predicted
Xp(:,:,1)=mvnrnd(proc_f(Xp(:,:,1))',proc_Q)';
%Posterior means of the particle filter
xfp(:,1)=Xp(:,:,1)*Wp(:,1);
%Posterior error covariance of the particle filter
Pfp(:,:,1)=(Xp(:,:,1)-xfp(:,1))*((Xp(:,:,1)-xfp(:,1))'.*Wp(:,1));
figure();
if bResample == true %resampling included
    Xr(:,:,1)=Xp(:,:,1); Wr(:,1)=Wp(:,1);
    for k=2:K+1
        [Xp(:,:,k),Wp(:,k)]=pfFilterStep(Xr(:,:,k-1),Wr(:,k-1)',Y(:,k-1),proc_f,proc_Q,meas_h,meas_R);
        if nargout > 2 
        xfp(:,k)=Xp(:,:,k)*Wp(:,k);
        Pfp(:,:,k)=(Xp(:,:,k)-xfp(:,k))*((Xp(:,:,k)-xfp(:,k))'.*Wp(:,k));
        end
        [Xr(:,:,k),Wr(:,k),j]=resampl(Xp(:,:,k),Wp(:,k)');
        plotFunc(k,Xp(:,:,k),Xp(:,:,k-1),j);
    end 
else %resampling decluded
for k=2:K+1
    [Xp(:,:,k), Wp(:,k)] = pfFilterStep(Xp(:,:,k-1), Wp(:,k-1)', ...
    Y(:,k-1), proc_f, proc_Q, meas_h, meas_R);
    if nargout>2 
        xfp(:,k)=Xp(:,:,k)*Wp(:,k);
        Pfp(:,:,k)=(Xp(:,:,k)-xfp(:,k))*((Xp(:,:,k)-xfp(:,k))'.*Wp(:,k));
    end
    plotFunc(k,Xp(:,:,k),Xp(:,:,k-1),[1:1:N:N]);
end
end
Xp=Xp(:,:,2:end); Wp=Wp(:,2:end);
if nargout > 2 
    xfp=xfp(:,2:end);
    Pfp=Pfp(:,:,2:end);
end
end
function [Xk, Wk, j] = resampl(Xk, Wk)
% Copy your code from previous task! 
[n,N]=size(Xk);
[Xr,j]=datasample(Xk,N,2,'Weights',Wk');
Wr=ones(1,N)./N;
end
function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
% Copy your code from previous task!
X_k=mvnrnd(proc_f(X_kmin1)',proc_Q)';
W_k=W_kmin1.*mvnpdf(yk',meas_h(X_k)',meas_R)';
W_k=W_k./sum(W_k);
end