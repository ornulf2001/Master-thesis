clc,clear

addpath(genpath("NewModel_2D"))
addpath(genpath("simulator"))

%Dynamics
run("parameters.m")
% z_eq=0.0365;
% Xeq = [0; z_eq; 0; 0; 0; 0];

index = @(A, i) A(i);
fz = @(z) index(f([0, z, zeros(1, 4)]', [0, 0]', params), 5);
zeq = fzero(fz, 0.1);
Xeq = [0, zeq, zeros(1, 4)]';
Ueq = [0, 0]';
[Ac, Bc, C] = linearizeModel(@f, @h, Xeq, Ueq, params);

%
clear params
run("parameters2.m")
delta = 1e-6;
Xeq = [0, 0.0365, zeros(1, 4)]';
Ueq = [0, 0]';
f2 = @(x,u) maglevSystemDynamics2d(x,u,params);
h2 = @(x,u) maglevSystemMeasurements2d(x,u,params);

[Ac2, Bc2, ~, ~] = finiteDifferenceLinearization(f2, h2, Xeq, Ueq, delta);

%Bc=-Bc;
C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
X0=[0;0;0;0;0;0];
NT=2000;

N_MHE=14;
N_MPC=50;
dt=0.0001;

R_MHE=1200*dt*diag([0.1,0.3,1]);
Q_MHE=1*diag([1,1,1,1,1,1]);
M_MHE=1*diag([3,3,3,1,1,1]);
noise_cov=0;

Q_MPC = diag([50000 1000 50000 0.5 0.1 0.5]);
R_MPC = diag([0.00005 0.00005]);

%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, nStates, nControls,MPC_options);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT);
MPC_Xopt2 = zeros(nStates, NT);
MPC_Uopt2 = zeros(nControls, NT);
MHE_est = zeros(nStates, NT);
yCurrent=zeros(nMeasurements,NT+1);  
yCurrent(:,1)= C*X0;
newY=yCurrent(:,1);
newU=[0;0];
MPC_Xopt(:, 1) = X0;
xCurrent = X0;
error=[];
for k=1:NT-1
    k
    
    
    if k<50
        U=[0;0];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xCurrent = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    elseif k>=50 && k<60
        U=[-10;-30];
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
    noise=noise_cov*dt*rand(size(C*xCurrent));
    yCurrent(:,k+1) = C*xCurrent + noise;
    newY=yCurrent(:,k+1);
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k)=xEst;
    error=[error,xCurrent-xEst];

end



% mhe=mhe.reset(X0);
% while ~mhe.isReadyToRun
%    newY=yCurrent(:,mhe.yBufferCount);
%    newU=MPC_Uopt(:,mhe.uBufferCount);
%    mhe=mhe.bufferInitialData(newY, newU);
% end
% xsol2=[];
% xsol2(:,1)=X0;
% for k=2:size(yCurrent,2)-1
%     newY=yCurrent(:,k);
%     newU=MPC_Uopt(:,k-1);
%     mhe=mhe.runMHE(newY,newU);
%     mhe.xCurrent
%     xsol2(:,k)=mhe.xCurrent;
% end



%Plot
figure(4)
plot(MHE_est(1,1:end)); hold on
plot(MPC_Xopt(1, 1:end));
legend('est', 'xopt')

figure(5)
plot(MHE_est(2,1:end)); hold on
plot(MPC_Xopt(2, 1:end));
legend('est', 'xopt')

figure(6)
plot(MHE_est(3,1:end)); hold on
plot(MPC_Xopt(3, 1:end));
legend('est', 'xopt')
% mpc.plotResults(MPC_Xopt,MPC_Uopt)














function [A,B,C] = linearizeModel(f,h,xlp,ulp,params)
    delta = 1e-10;
    
    A = zeros(6, 6);    
    for i = 1:6
        xp = xlp;
        xn = xlp;
    
        xp(i) = xp(i) + delta;
        xn(i) = xn(i) - delta;
        
        dxp = f(xp, ulp, params);
        dxn = f(xn, ulp, params);
        
        A(:, i) = (dxp - dxn)/(2*delta);
    end
    
    B = zeros(6, 2);    
    for i = 1:2
        up = ulp;
        un = ulp;
    
        up(i) = up(i) + delta;
        un(i) = un(i) - delta;
        
        dxp = f(xlp, up, params);
        dxn = f(xlp, un, params);
        
        B(:, i) = (dxp - dxn)/(2*delta);
    end
    
    C = zeros(2, 6);    
    for i = 1:6
        xp = xlp;
        xn = xlp;
    
        xp(i) = xp(i) + delta;
        xn(i) = xn(i) - delta;
        
        yp = h(xp, params);
        yn = h(xn, params);
        
        C(:, i) = (yp - yn)/(2*delta);
    end
end
