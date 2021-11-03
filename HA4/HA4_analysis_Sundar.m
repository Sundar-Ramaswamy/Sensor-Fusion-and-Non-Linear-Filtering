%% 1A
clear all; 
close all; 
clc;
cp = fp.getColor(1:10);

sigma_v = 1*1e-4;
sigma_w = pi/180;
%% True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
% Initial state
x_0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x_0),K+1);
X(:,1) = x_0;
% Create true track
for i=2:K+1
% Simulate
X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
% Set turnrate
X(5,i) = omega(i);
end

% Prior information
x_0 = [0 0 0 0 0]';
P_0 = diag([10 10 10 5*pi/180 pi/180].^2);
% Sensor positions
S1 = [300 -100]';
S2 = [300 -300]';

% measurement noise standard deviation
R = diag([pi/180 pi/180].^2);
% generate measurement sequence
h = @(x) rangeBearingMeasurements(x, S1);
Y = genNonLinearMeasurementSequence(X,h,R);
% Motion model
f = @(x) coordinatedTurnMotion(x,T);
Q = diag([0 0 T*sigma_v^2 0 T*sigma_w^2]);

%[xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, @sigmaPoints, 'CKF');

X_m = S1 + Y(1,:).*[cos(Y(2,:));sin(Y(2,:))];
figure
subplot(1,2,1)
title('Filter')
plotTurnU(X, xf, Pf, X_m, S1,'filter',2)
subplot(1,2,2)
title('Smoother')
plotTurnU(X, xs, Ps, X_m, S1,'smoother',4)
figure
plotTurnUError( T, xf, xs, X )
function plotTurnU( X, xf, Pf, Xm, s1, signame, coln)
    cp = fp.getColor(1:10);
    grid on; hold on, axis equal;
    for i=1:5:length(xf)
        ell_xy = sigmaEllipse2D(xf(1:2,i),Pf(1:2,1:2,i),3,50);
        p5 = fill(ell_xy(1,:),ell_xy(2,:), cp(coln,:),'facealpha',.1, 'DisplayName',[signame,' 3-sigma']);   %,'edgecolor','none'
    end

    p1 = plot(X(1,:),X(2,:),   '-', 'Color', cp(1,:), 'LineWidth',2, 'DisplayName','True position sequence');
    p2 = plot(xf(1,:),xf(2,:), '-', 'Color', cp(coln,:), 'LineWidth',2, 'DisplayName',[signame,' position']);
    sc1 = scatter(s1(1), s1(2), 100, 'o', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor', cp(4,:), 'MarkerEdgeColor', cp(4,:),'DisplayName','sensor 1 location');

    axis manual
    p4 = plot(Xm(1,:),Xm(2,:), 'Color', [cp(3,:) 0.3], 'LineWidth',1, 'DisplayName','Measured position');

    xlabel 'pos x', ylabel 'pos y'
    legend([p1 p2 p4 sc1 p5], 'Location','west')
    % if savefig fp.savefig(sprintf('q3_%s',name2save)); end
end

function plotTurnUError( T, xf, xs, X )
    cp = fp.getColor(1:10);
    K = length(X)-1;
    grid on, hold on;
    p1 = plot( (1:K)*T, vecnorm(xf(1:2,:)-X(1:2,2:end), 2, 1), 'Color', cp(2,:) , 'LineWidth',2, 'DisplayName','Filter error');
    p2 = plot( (1:K)*T, vecnorm(xs(1:2,:)-X(1:2,2:end), 2, 1), 'Color', cp(4,:) , 'LineWidth',2, 'DisplayName','Smoother error');
    ylabel('$|p_k - \hat{p_{k|k}}|_2$', 'Interpreter','Latex', 'FontSize',16), xlabel('Time [s]')
    title 'Position error'
    legend([p1 p2])
    % if savefig fp.savefig(sprintf('q3_%s_err',name2save)); end
end
