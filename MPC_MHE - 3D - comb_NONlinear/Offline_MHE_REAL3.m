%%
clc,clear

addpath(genpath('3D model reduced order_fixed'))
%addpath(genpath("datasets"))

load("data_dMag2.mat")

%Dynamics
params = parameters1();

index = @(A, i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);   
zeq =  fzero(fz,0.03);
xeq = [0,0,zeq,zeros(1,7)]';
xlp=xeq;
ueq = [0,0,0,0]';
ulp=ueq;
X0=xlp;
%[Ac, Bc, C] = linearizeModel(@f, @h, xeq, ueq, params);
load("ABC_simple_model_reduced.mat");Ac=A; Bc=B;

%% Load data
I = 100:1400;
u = [data.u.Ix_plus(I), data.u.Iy_plus(I), data.u.Ix_minus(I), data.u.Iy_minus(I)]';
y = 1e-3*[
    data.sensorData{1}.bx(I), data.sensorData{1}.by(I), data.sensorData{1}.bz(I),...
    data.sensorData{2}.bx(I), data.sensorData{2}.by(I), data.sensorData{2}.bz(I),...
    data.sensorData{3}.bx(I), data.sensorData{3}.by(I), data.sensorData{3}.bz(I)
]';
yeq    = mean(y(:,end-50:end),2);
t      = data.t(I);
t      = t-t(1);
dtvec  = diff(data.t);
dt     = dtvec(1);
%% MHE Tuning
alpha=1;
N_MHE=10;



%noise_std=0.1*1e-3; %mT
%R_MHE=inv(noise_std^2*eye(size(C,1)));  %Measurement noise weight = inv(measurement noise cov)  

%Q_MHE=1e8*eye(size(Ac,1));
%Q_MHE=1e6*diag([1e6,1e5,1e5,1e5,1e5,0.5e-3,0.5e-3,0.5e-3,0.5e-3,0.5e-3]);  
Q_MHE=1e6*diag([1e1,1e1,1e1,1e1,1e1,1e5,1e5,1e5,1e5,1e5]);  
R_MHE=inv(cov(y(:,400:end)'));

M_MHE = 1e2*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]); %Arrival cost weight initial guess (updates KF-style in loop)
%M_MHE = 1e2*eye(size(Ac,1));
P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling =1e-4; %Scaling the weight matrices uniformly to ease solving

%% Run
MHE_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);
NT=ceil(size(y,2));

% Initialization
xhat=X0;
Phat=P0;
KF_est=zeros(size(Ac,1),NT-1);
KF_est(:,1)=xhat;
MHE_est=zeros(size(Ac,1),NT-1);
vsol=zeros(mhe.nMeasurements,NT-1);
wsol=zeros(size(Ac,1),NT-1);
MHE_est(:,1)=X0-xlp;
newY_f=y(:,1);

% KF Tuning
R_KF= inv(R_MHE);
%Q_KF = 1e-8*eye(size(Ac));
%Q_KF= inv(1e6*diag([1e4,1e4,1e4,1e4,1e4,0.5e-3,0.5e-3,0.5e-3,0.5e-3,0.5e-3]));
Q_KF=inv(Q_MHE);

A = expm(Ac * dt);
B = (A - eye(size(Ac))) * (Ac \ Bc);

for k=1:NT-1
    k

    newY     = y(:,k+1);
    newY_f   = alpha*newY + (1-alpha)*newY_f; %EMA prefilter before MHE
    newY     = newY_f;
    newU     = u(:,k);

    mhe           = mhe.runMHE(newY- yeq, newU - ueq);
    MHE_est(:,k+1)  = mhe.xCurrent;
    vsol(:,k+1)   = mhe.vCurrent;
    wsol(:,k)     = mhe.wCurrent;
    

    xhat_pred  = A*xhat + B*(newU-ueq);
    Phat_pred  = A*Phat*A' + Q_KF;

    y_pred     = C*xhat_pred;
    y_unbiased = y(:,k+1) - yeq;
    innovation = y_unbiased - y_pred;

    S = C*Phat_pred*C' + R_KF; 
    K = Phat_pred*C'/S;    

    xhat = xhat_pred + K*innovation;
    Phat = (eye(size(Phat_pred)) - K*C)*Phat_pred*(eye(size(Phat_pred)) - K*C)' + K*R_KF*K'; %Joseph form 
    
    KF_est(:, k+1) = xhat;
 
end


%% Reconstructing absolute state estsimates and estimated outputs 
MHE_est = MHE_est + xlp;
KF_est  = KF_est + xlp;

y_est_KF  =  yeq + C*KF_est ; % KF
y_est_MHE =  C*MHE_est+yeq; % MHE

%% Comparing dy_est vs IIR dy
dyx0  = 1e-3* data.dMagField.dMagFieldX(I)';
dyy0  = 1e-3* data.dMagField.dMagFieldY(I)';

for i=1:NT
    dxEst(:,i) = circshift(MHE_est(:,i), 5);
    dyEst(:,i) = C*dxEst(:,i);
end

x_recon = zeros(size(MHE_est,1)/2,size(MHE_est,2));  % Reconstructed x from dxEst
for j=1:size(Ac,1)/2
    x_recon(j,:) = MHE_est(j,60) + cumtrapz(t, MHE_est(j+5,:));
end

diffY=diff(y(1,:))./diff(t)';

figure(9);clf
plot(dyx0,"b-"); hold on
plot(dyEst(1,2:end),"r-"); 
plot(diffY,"k-")
legend(["dMagFieldX","dY\_X\_est"])
legend(["dMagFieldX","filter","meas"])

figure(10);clf
plot(t,x_recon(1:3,:));hold on
title("MHE estimates position reconstructed from dxEst")
yline(zeq,"r--")
yline(0,"k--")
legend(["x","y","z","zeq"])


%% Plotting
figure(1)
clf
h1 = plot(t/dt, y(1:6,:), 'b-'); hold on
h2 = plot(t/dt, y_est_KF(1:6,:), 'r--');
legend([h1(1), h2(1)], ["meas", "est"]) 
title("KF")
ylim([1.25*min(y(:)), 1.25*max(y(:))])

figure(2)
clf
h1 = plot(t/dt, y(1:6,:), 'b-'); hold on
h2 = plot(t/dt, y_est_MHE(1:6,:), 'r--');
legend([h1(1), h2(1)], ["meas", "est"]) 
title("MHE")
ylim([1.25*min(y(:)), 1.25*max(y(:))])

figure(3);clf
plot(t,MHE_est(1:3,:));hold on
title("MHE estimates position")
yline(zeq,"r--")
yline(0,"k--")
legend(["x","y","z","zeq"])

figure(4);clf
plot(t,MHE_est(6:8,:));hold on
title("MHE estimates velocity")
yline(0,"k--")
legend(["xdot","ydot","zdot"])

figure(5);clf
plot(t,KF_est(1:3,:));hold on
title("KF estimates position")
yline(zeq,"r--")
yline(0,"k--")
legend(["x","y","z","zeq"])

figure(6);clf
plot(t,KF_est(6:8,:));hold on
title("KF estimates velocity")
yline(0,"k--")
legend(["xdot","ydot","zdot"])

figure(7);clf
plot(t(1:end-1),wsol(:,:))
yline(0,"k--")
legend(["x","y","z","ph","th","xdot","ydot","zdot","phdot","thdot"])
title("MHE process noise w estimates")

figure(8);clf
plot(t,vsol(:,:))
yline(0,"k--")
legend(["x","y","z","ph","th","xdot","ydot","zdot","phdot","thdot"])
title("MHE measurement noise v estimates")



%%
function data = getParams()
    persistent params

    if isempty(params)
        parameters2;
    end
    data = params;
end