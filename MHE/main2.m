clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*

%Dynamics
params=getParams();
z_eq=0.0365;
Xeq = [0; z_eq; 0; 0; 0; 0];
Ueq=[0;0];
delta = 1e-6;
f = @(x,u) maglevSystemDynamics2d(x,u,params);
h = @(x,u) maglevSystemMeasurements2d(x,u,params);
[Ac, Bc, ~, ~] = finiteDifferenceLinearization(f, h, Xeq, Ueq, delta);
C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
X0=[0;0;0;0;0;0];
NT=400;

N_MHE=3;
N_MPC=10;
dt=0.003;

R_MHE=5*dt*diag([1,0.1,1]);
Q_MHE=diag([1,1,1,1,1,1]);
M_MHE=0.1*diag([1,1,1,3,3,3]);
noise_cov=0.05;

Q_MPC = diag([50 100 50 0.8 0.5 0.8]);
R_MPC = diag([0.01 0.1]);

%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, nStates, nControls,MPC_options);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
yNext=zeros(nMeasurements,NT+1);  
yNext(:,1)= C*X0;
newY=yNext(:,1);
xNext = X0;
MPC_Xopt(:, 1) = X0;
error=[];
for k=1:NT-1
  
    if k<50
        U=[0;0];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xNext = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    elseif k>=50 && k<60
        U=[-10;-30];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xNext = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    else
        
        [xNext, Uopt]=mpc.runMPC(xEst);
        MPC_Xopt(:,k+1) = xNext;
        MPC_Uopt(:,k)= Uopt;
        newU=MPC_Uopt(:,k);
    end
    noise=noise_cov*dt*randn(size(C*xNext));
    yNext(:,k+1) = C*xNext + noise;
    newY=yNext(:,k+1)
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k)=xEst;
    error=[error,xNext-xEst];

end




% figure(4)
% plot(MHE_est(1,1:end)); hold on
% plot(MPC_Xopt(1, 1:end));
% plot(yNext(1,:))
% title("X")
% legend('est', 'xopt',"meas")
% 
% figure(5)
% plot(MHE_est(2,1:end)); hold on
% plot(MPC_Xopt(2, 1:end));
% title("Z")
% legend('est', 'xopt')
% 
% figure(6)
% plot(MHE_est(3,1:end)); hold on
% plot(MPC_Xopt(3, 1:end));
% title("Theta")
% legend('est', 'xopt')
mpc.plotResults(MPC_Xopt,MPC_Uopt)