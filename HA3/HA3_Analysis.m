close all; 
clear all; 
clc;
%% Problem 1: Approximations of mean and covariance
% number of samples
N=10000;
%update type
type = 'CKF';
% two prior state_densitiess 1 or 2
state_densities = '1';  
switch state_densities
case '1'
X_0_mean = [125 125]';
P_0 = [100 0;0 25];
case '2'
X_0_mean = [-25 125]';
P_0 = [100 0;0 25];
case '3'
X_0_mean = [60 60]';
P_0 = [100 0;0 25];
end
%sensor position
S1 = [0;100]';
S2 = [100;0]';
% Noise covariance
R = diag([0.1*pi/180 0.1*pi/180].^2);
% Measurement model from functions
h = @(x) dualBearingMeasurement(x,S1,S2);
MeasurementSequence = @(x)genNonLinearMeasurementSequence(x, h, R);
%1a
[y_mu, y_sigma, y_s] = approxGaussianTransform(X_0_mean, P_0, @(x)genNonLinearMeasurementSequence(x,h,R) , N);
%1c
switch type
case 'UKF'
[SP,W] = sigmaPoints(X_0_mean,P_0,type);
hSP1 = h(SP);
[ye_mu, ye_sigma] = nonLinKFprediction(X_0_mean, P_0, h, R, type);
case 'CKF'
[SP,W] = sigmaPoints(X_0_mean,P_0,type);
hSP1 = h(SP);
[ye_mu, ye_sigma] = nonLinKFprediction(X_0_mean, P_0, h, R, type);
case 'EKF'
[hx, dhx] = h(X_0_mean);
[ye_mu, ye_sigma] = nonLinKFprediction(X_0_mean, P_0, h, R, type);
end

level=3;
[xy1] = sigmaEllipse2D(y_mu, y_sigma, level,200);
[xy2] = sigmaEllipse2D(ye_mu, ye_sigma, level,200);
figure(1);
scatter(y_s(1,:), y_s(2,:),5,'*r')
hold on
plot(xy1(1,:), xy1(2,:),'Linewidth',2);
scatter(y_mu(1), y_mu(2),25, 'og');
switch type
case 'UKF'
scatter(hSP1(1,:), hSP1(2,:),'om');
case 'CKF'
scatter(hSP1(1,:), hSP1(2,:),'ok');
end
hold on
scatter(ye_mu(1,:), ye_mu(2,:),'oc')
plot(xy2(1,:),xy2(2,:))
xlabel('phi_1');
ylabel('phi_2');
switch type
case 'UKF'
legend('Measurement state densities','Untransformed ellipse','untransformed mean','sigma points','transformed mean','transformed ellipse','Location','best')
case 'CKF'
legend('Measurement state densities','Untransformed ellipse','untransformed mean','sigma points','transformed mean','transformed ellipse','Location','best')
case 'EKF'
legend('Measurement state densities','Untransformed ellipse','untransformed mean','transformed mean','transformed ellipse','Location','best')
end

%% problem 2: Non-linear Kalman 
x_0=[0 0 20 0 (5*pi)/180]';
P_0= diag([100 100 4 (pi/180)^2 (pi/180)^2]);
s1=[-200,100]';
s2=[-200,-100]';
T=1;
N=100;
sigma_v = 1;
sigma_w = (pi/180);
if state_densities == 1
    sigma_phi_1 = (10*pi/180);
else
    sigma_phi_1 = (0.5*pi/180);
end
sigma_phi_2 = (0.5*pi/180);
Q = diag([0, 0, T*sigma_v, 0, T*sigma_w].^2);
R = diag([sigma_phi_1, sigma_phi_2].^2);
% Process and Measurement Model and Transformed Sequence
f=@(x)coordinatedTurnMotion(x,T);
h=@(x)dualBearingMeasurement(x,s1,s2);
X=genNonLinearStateSequence(x_0, P_0, f, Q, N);
Y=genNonLinearMeasurementSequence(X, h, R);
Xm(1,:) = (s2(2)-s1(2)+tan(Y(1,:))*s1(1)- tan(Y(2,:))*s2(1))./(tan(Y(1,:))- tan(Y(2,:)));
Xm(2,:) = s1(2)+tan(Y(1,:)).*(Xm(1,:)- s1(1));

%2a and 2b
level=3;
for type = {'EKF','UKF','CKF'}
    [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type{1});  
    figure();
    clf;
    hold on
    plot(s1(1),s1(2),'*r','Linewidth',2);
    plot(s2(1),s2(2),'or','Linewidth',2);
    plot(X(1,:),X(2,:),'-m');
    plot(xf(1,:),xf(2,:),'-.c');
    plot(Xm(1,:),Xm(2,:),'+g');
	for i= 1:5:length(xf)
        [xy3]=sigmaEllipse2D(xf(1:2,i),Pf(1:2,1:2,i),level,50);
        plot(xy3(1,:),xy3(2,:),'--b');
	end
    xlabel('x');
    ylabel('y');
    legend('Sensor 1 Pos','Sensor 2 Pos','True State','Filtered State','Measured State','3sigma region','Interpreter','Latex','Location','best');
    hold off     
end

% Problem 2.c: Plot Histograms of Estimation Error
MC = 100;
est_err = cell(1,3);
type = {'EKF','UKF','CKF'};
for imc = 1:MC    
    X = genNonLinearStateSequence(x_0, P_0, f, Q, N);
    Y = genNonLinearMeasurementSequence(X, h, R);
     for itype = 1:numel(type)
        
        %  Kalman filter
        [xfM,PfM,xpM,PpM] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,type{itype});
        est_err{1,itype}(1:2,end+1:end+length(xf)) = X(1:2,2:end) - xfM(1:2,:);

    end
end

 MCcount = 100;
 close all;
 loc = {'x','y'};
     figure(2);
     for itype = 1:numel(type)
        for iloc = 1:numel(loc)
            subplot(2,3, itype + (iloc-1)*numel(type) );
            hold on;
%             
             histo = est_err{1,itype}(iloc,:);
             errmu  = mean(histo);
             errsig = std(histo);

            histogram( histo, MCcount ,'Normalization','pdf');
            level=3;
            N2=100;
            x = linspace(errmu-level*sqrt(errsig^2), errmu+level*sqrt(errsig^2), N2);
            y = normpdf(x, errmu, sqrt(errsig^2));
            plot(x,y, 'LineWidth',2 );
   
        end
     end
        
%% Problem 3: Tuning non-linear filters
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
% Simulate
X(:,i) = coordinatedTurnMotion(X(:,i-1),T);
% Set turn-rate
X(5,i) = omega(i);
end
% Prior
x_0 = [0 0 0 0 0]';
P_0 = diag([10 10 10 5*pi/180 pi/180].^2);
% Sensor positions
s_1 =[300 -100]';
s_2 =[300 -300]';
gamma=[0 0;0 0;1 0;0 0;0 1];
Q=gamma*diag([200*1 200*pi/180].^2)*gamma';
% measurement variance
R = 0.05*diag([pi/180 pi/180].^2);
% generate measurement sequence
h = @(x) dualBearingMeasurement(x,s_1,s_2);
Y = genNonLinearMeasurementSequence(X,h,R);
% Motion model
motionModel=@(x)coordinatedTurnMotion(x,T);
[xf,Pf,xp,Pp]=nonLinearKalmanFilter(Y,x_0,P_0,motionModel,Q,h,R,'CKF');
%unfiltered position
xmeas =(s_2(2)-s_1(2)+tan(Y(1,:))*s_1(1)-tan(Y(2,:))*s_2(1))./(tan(Y(1,:))-tan(Y(2,:)));
ymeas =s_1(2)+tan(Y(1,:)).*(xmeas(1,:) - s_1(1));
figure(1);
grid on;
hold on
axis equal;
plot(X(1,:),X(2,:),'m');
plot(xf(1,:),xf(2,:),'r');
scatter(s_1(1),s_1(2),100,'o');
scatter(s_2(1),s_2(2),200,'o');
axis manual
plot(xmeas,ymeas,'*');
for i=1:15:length(xf)
variance_xy = sigmaEllipse2D(xf(1:2,i),Pf(1:2,1:2,i),3,50);
plot(variance_xy(1,:),variance_xy(2,:))
end
xlabel('pos x');
ylabel('pos y');
legend('true state','filtered position','sensor1 potion','sensor2 potion','Measurements')
% plot position error
err=X(1:2,2:end) - xf(1:2,:);
figure(2);
grid on;
hold on;
plot((1:K)*T,err(1,:),(1:K)*T,err(2,:))

