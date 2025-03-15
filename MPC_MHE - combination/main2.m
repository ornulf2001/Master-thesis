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
X0=[0.001;0.001;0;0;0;0];
NT=500;

N_MHE=10;
N_MPC=20;
dt=0.003;

alpha=0.7;
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));         
Q_MHE=2e10*diag([10,10,10,1,1,1]);
M_MHE=5e5*diag([10,10,10,3,3,3]);
P0 = inv(M_MHE);

Q_MPC = diag([100 100 100 3 3 3]);
R_MPC = diag([0.001 0.001]);

Q_LQR = diag([10,10,10,1,1,1]);
R_LQR = diag([0.01,0.01]);
%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,1e-5*Q_MHE,1e-5*R_MHE,1e-5*M_MHE,X0,P0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, nStates, nControls,MPC_options, xRef, lbuRef, ubuRef);

MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
MHE_est(:,1)=mhe.x0; xEst = mhe.x0;
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
    if k<=60
        [K_dlqr,~,~] = dlqr(mpc.A, mpc.B, Q_LQR, R_LQR);
        K = [120,120,120,5,5,5;120,120,120,5,5,5];
            %Kp = K_dlqr(:,1:nStates/2);
            %Kd = K_dlqr(:,nStates/2+1:nStates);
            %U_PD = Kp*(xEst(1:nStates/2)-xRef(1:nStates/2)) + Kd*(xEst(nStates/2+1:nStates)-xRef(nStates/2+1:nStates));
        U_LQR = -K_dlqr*(xEst-xRef);

        U=U_LQR;

        fade_period=50:60;
        if ismember(k,fade_period)
            [~, Uopt]=mpc.runMPC(xEst);
            Udiff=U_LQR - Uopt;
            U = (1-gamma_f(k,fade_period))*U_LQR + gamma_f(k,fade_period)*Uopt;
        end
        
        %[T, X] = ode15s(@(t, x) f(x, U, params), tspan, MPC_Xopt(:,k));
        %MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;

        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    else
        [~, Uopt]=mpc.runMPC(xEst);

        if k>=65 && k<75    
            Uopt=Uopt+[50;150];
        end

        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*(Uopt+uRef);
        %[T, X] = ode15s(@(t, x) f(x, Uopt, params), tspan, MPC_Xopt(:,k));
        %MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        MPC_Uopt(:,k)= Uopt + uRef;
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
%mpc.plotResults(MPC_Xopt,MPC_Uopt,xRef)

% figure(7)
% plot(yNext(1,:)); hold on
% plot(yNext_f(1,:))


function gamma = gamma_f(k,fade_period)
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end