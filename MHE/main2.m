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
NT=600;

N_MHE=8;
N_MPC=30;
dt=0.002;

R_MHE=1*dt*diag([0.1,0.1,0.1]);
Q_MHE=diag([1,1,1,1,1,1]);
M_MHE=0.2*diag([1,1,1,3,3,3]);
noise_cov=0.05;

Q_MPC = diag([50 300 50 0.1 0.5 0.3]);
QN_MPC =1000*Q_MPC;
R_MPC = diag([0.05 0.001]);

%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, QN_MPC, nStates, nControls,MPC_options);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT);
MHE_est = zeros(nStates, NT);
yCurrent=zeros(nMeasurements,NT+1);  
yCurrent(:,1)= C*X0;
newY=yCurrent(:,1);
newU=[0;0];
MPC_Xopt(:, 1) = X0;
xCurrent = X0;

for k=1:NT
    tic
    k
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k)=mhe.xCurrent;
    if k<20
        U=[0;0];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xCurrent = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    elseif k>=20 && k<22
        U=[-5;-13];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xCurrent = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);
    else
        [xNext, Uopt]=mpc.runMPC(xEst);
        MPC_Xopt(:,k+1) = xNext;
        MPC_Uopt(:,k)= Uopt;
        xCurrent=xNext;
        newU=MPC_Uopt(:,k);
    end

    yCurrent(:,k+1) = C*xCurrent + noise_cov*dt*rand(size(C*xCurrent));
    newY=yCurrent(:,k+1);
    toc
end



mpc.plotResults(MPC_Xopt,MPC_Uopt)