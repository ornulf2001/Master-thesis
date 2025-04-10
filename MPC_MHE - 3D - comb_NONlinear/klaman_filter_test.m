%% Load data
load('data_with_control_fixed_double_sample_rate.mat');

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
y_eq = h(xeq, ueq, params);

%% Kalman Filter
% Dimensions
n = length(xeq);
m = size(y, 2);
N = length(t);

% Initialization
xhat = zeros(n, N);
xhat(:,1) = xeq;

P    = zeros(n, n, N);
P(:,:,1)  = 1e-2*eye(n);

% Covariances
Q = 1e-2*eye(n);
R = diag(diag(cov(y(400:end,:))));

% Loop
for k = 1:N-1
    dt = t(k+1) - t(k);
    
    % Discretization
    A_d = eye(n) + A*dt;
    B_d = B*dt;
    
    % Prediction Step:
    xhat_minus = A_d*xhat(:,k) + B_d*(u(k,:)' - ueq);
    P_minus = A_d*P(:,:,k)*A_d' + Q;
    
    % Measurement Prediction (Adding yeq instead of subtracting)
    ypred = y_eq + C*(xhat_minus - xeq);
    
    % Innovation
    innovation = y(k+1,:)' - ypred;
    
    % Kalman Gain computation
    K = P_minus*C'/(C*P_minus*C' + R);
    
    % Update Step
    xhat(:,k+1) = xhat_minus + K*innovation;
    P(:,:,k+1)  = (eye(n) - K*C)*P_minus;
end

% Recompute estimated y
yhat = repmat(y_eq, 1, N) + C*(xhat - repmat(xeq, 1, N));

%% Plot measured vs estimated outputs
figure(2);
clf; hold on; box on;
plot(t,y(:,1:6), 'b-')
YLIM = ylim();
plot(t,yhat(1:6,:)', 'r--')
ylim(YLIM)
xlabel('Time (s)')
ylabel('Sensor Measurements')
legend({'Measured y','','','Estimated yhat'})

figure(3);
clf; 
subplot(2,1,1);
hold on; box on;
plot(t,xhat(1:3,:)', 'b')
xlabel('Time (s)')
ylabel('Estimated positions')

subplot(2,1,2);
hold on; box on;
plot(t,xhat(4:6,:)', 'b')
xlabel('Time (s)')
ylabel('Estimated velocities')
