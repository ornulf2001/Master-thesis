clc,clear

addpath(genpath("NewModel_2D"))
addpath(genpath("simulator"))

%Dynamics
run("parameters.m")

index = @(A, i) A(i);
fz = @(z) index(f([0, z, zeros(1, 4)]', [0, 0]', params), 5);
zeq = fzero(fz, 0.1);
Xeq = [0, zeq, zeros(1, 4)]';
Ueq = [0, 0]';
[Ac, Bc, C] = linearizeModel(@f, @h, Xeq, Ueq, params);

%C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
xRef = [0; 0.03; 0; 0; 0; 0];
X0=[0;0;0;0;0;0];
NT=500;

N_MHE=15;
N_MPC=20;
dt=0.003;

alpha=0.8;
noise_std=0.1*1e-3; %mT
R_MHE=1e-5*inv(noise_std^2*eye(nMeasurements));         %5*dt*diag([1,1]);
Q_MHE=1e-5*2e8*diag([10,1,10,10,1,10]);
M_MHE=1e-5*5e6*diag([10,1,10,3,1,3]);
P0 = inv(M_MHE);

Q_MPC = diag([100 50 100 1 1 1]);
R_MPC = diag([0.003 0.001]);
%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, nStates, nControls,MPC_options, xRef, lbuRef, ubuRef);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
MHE_est(:,1)=mhe.x0;
yNext=zeros(nMeasurements,NT+1);  
yNext(:,1)= C*X0;
yNext_f=zeros(nMeasurements,NT+1);
yNext_f(:,1)=C*X0;
newY=yNext(:,1);
xNext = X0;
MPC_Xopt(:, 1) = X0;
error=[];

tspan = [0, dt];

% Calculating the reference input for stabilizing in the reference point
uRef = mpc.computeReferenceInput();

%init_plot(X0);
for k=1:NT-1
    k
    if k<50
        U=[0;0];
        %MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        [T, X] = ode15s(@(t, x) f(x, U, params), tspan, MPC_Xopt(:, k));
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        %xNext = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    else
        xEst
        MPC_Xopt(:,k)
       [xNext, Uopt]=mpc.runMPC(xEst);
        %MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*Uopt;
        [T, X] = ode15s(@(t, x) f(x, Uopt, params), tspan, xEst);
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        MPC_Uopt(:,k)= Uopt + uRef;
        newU=MPC_Uopt(:,k);
    end
        
    if k>=60 && k<70    

        Uopt=Uopt+[10;30];
        %MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*Uopt;
        [T, X] = ode15s(@(t, x) f(x, U, params), tspan, xEst);
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        %xNext = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = Uopt;
        newU=MPC_Uopt(:,k);
    end
    
    noise=noise_std*randn([nMeasurements,1]);
    yNext(:,k+1) = C*MPC_Xopt(:,k+1)+ noise;
    yNext_f(:,k+1)=alpha*yNext(:,k+1) + (1-alpha)*yNext_f(:,k);
    newY=yNext_f(:,k+1);
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k+1)=xEst;
    %error=[error,xNext-xEst];

end



%Plot
figure(4)
%plot(yNext_f(1,:)); hold on
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
%mpc.plotResults(MPC_Xopt,MPC_Uopt)

% figure(7)
% plot(yNext(1,:)); hold on
% plot(yNext_f(1,:))