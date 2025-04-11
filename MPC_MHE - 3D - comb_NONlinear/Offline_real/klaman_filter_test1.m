%% Load data
clc,clear
load('data_dMag2.mat');

%%
I = 100:700;

t = data.t(I);
t = t - t(1);
u = [data.u.Ix_plus(I), data.u.Iy_plus(I), data.u.Ix_minus(I), data.u.Iy_minus(I)];
y = 1e-3*[
    data.sensorData{1}.bx(I), data.sensorData{1}.by(I), data.sensorData{1}.bz(I),...
    data.sensorData{2}.bx(I), data.sensorData{2}.by(I), data.sensorData{2}.bz(I),...
    data.sensorData{3}.bx(I), data.sensorData{3}.by(I), data.sensorData{3}.bz(I)
];

%% Load model
params = parameters;

%% Find equilibrium point
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);

zeq = fzero(fz, 0.1);

xeq = [0,0,zeq,zeros(1,7)]';
ueq = [0,0,0,0]';

%% Linearize model
xlp = xeq;
ulp = ueq;

[A,B,C] = linearizeModel(@f, @h, xlp, ulp, params);

% Compute equilibrium measurement
yeq = h(xeq, ueq, params);
yeq = y(end,:)'; % HACK

%% Kalman Filter
% Dimensions
n = length(xeq);
m = size(y, 2);
N = length(t);

% Initialization
xhat = zeros(n, N);
xhat(:,1) = xeq;

P    = zeros(n, n, N);
P(:,:,1)  = 1e0*eye(n);

% Covariances
Q = 1e-5*eye(n);
R = diag(diag(cov(y(400:end,:))));

% Loop
for k = 1:N-1
    dt = t(k+1) - t(k);
    
    % Discretization
    sysd = c2d(ss(A,B,C,0), dt, 'zoh');
    A_d = sysd.A;
    B_d = sysd.B;

    % Prediction Step:
    xhat_minus = A_d*xhat(:,k) + B_d*(u(k,:)' - ueq);
    P_minus = A_d*P(:,:,k)*A_d' + Q;
    
    % Measurement Prediction
    ypred = C*xhat_minus;
    
    % Innovation
    y_unbiased = y(k+1,:)' - yeq;
    innovation = y_unbiased - ypred;
    
    % Kalman Gain computation
    K = P_minus*C'/(C*P_minus*C' + R);
    
    % Update Step
    xhat(:,k+1) = xhat_minus + K*innovation;
    P(:,:,k+1)  = (eye(n) - K*C)*P_minus;
end

% Recompute estimated y
yhat = C*xhat + yeq;
xhat = xhat + xeq;

%% Plot measured vs estimated outputs
figure(2);
clf; hold on; box on;

plot(t,y(:,1:6), 'b-')

%YLIM = ylim();

plot(t,yhat(1:6,:)', 'r--')
yline(yeq(1:6), 'g-', 'LineWidth', 2)

%ylim(YLIM)

xlabel('Time (s)')
ylabel('Sensor Measurements')
legend({'Measured y','','','Estimated yhat'})

figure(3);
clf; hold on; box on;
plot(t,xhat(1:3,:))

figure(4);
clf; hold on; box on;
plot(t,xhat(6:8,:))