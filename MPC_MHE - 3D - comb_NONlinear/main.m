clc,clear
%addpath(genpath("3D model"))
addpath(genpath('3D model reduced order'))


% Define system parameters
params = parameters;

% Find equilibrium point
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1

zeq =  fzero(fz,0.1)

xeq = [0,0,zeq,zeros(1,7)]';
ueq = [0,0,0,0]';

% Linearize model
xlp = xeq;   % 10x1 equilibrium state
ulp = ueq;

[Ac,Bc,C] = linearizeModel(@f,@h,xlp,ulp,params);


nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%% Tuning
xRef = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
X0=[0.001;0.001;zeq;0;0;0;0;0;0;0;];
NT=500;

N_MHE=16;
N_MPC=20;
dt=0.003;


% States: | x y z phi theta xdot ydot zdot phidot thetadot |

alpha=0.7;
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));         
Q_MHE=25e4*diag([100,100,100,100,100,100,100,100,100,100]); 
    %Start out with low Q to trust measurements during start up, 
    %then increase Q after N_MHE+1. 
    %See below in loop
                                  
%load("KF_M.mat")
%M_MHE = mheM;
% M_MHE=5e5*diag([100,100,1000,100,100,1000,1000,30,30,30]);
M_MHE = 1e5*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
P0 = inv(M_MHE);

Q_MPC = diag([50 50 800 10 10 1 1 100 1 1]);
R_MPC = diag([0.08 0.08 0.08 0.08]);

%Q_LQR=Q_MPC;
%R_LQR=R_MPC;

%Q_LQR = diag([500 500 10 0.8 0.8 7 0.8 0.8 0.8 0.8]);
%R_LQR = diag([0.001 0.001 0.0001 0.0001]);
Q_LQR = diag([ ...
   1e1,1e1,1e1,1e1,1e1, ...
   1e1,1e1,1e5,1e1,1e1
   ]);
R_LQR = 1e2*eye(4);


% Bounds
run("mpc_bounds.m")

% Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,1e-8*Q_MHE,1e-8*R_MHE,1e-8*M_MHE,X0,xlp,P0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, [], [], Q_MPC, R_MPC, nStates, nControls,MPC_options, xRef, [], []);

X_sim = zeros(nStates, NT);
U_sim = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
MHE_est(:,1)=mhe.x0; xEst = mhe.x0;
yNext=zeros(nMeasurements,NT);  
yNext(:,1)= C*(X0-xlp);
% yNext(:, 1) = h(X0, params);
yNext_f=zeros(nMeasurements,NT);
yNext_f(:,1)=C*(X0-xlp);
% yNext_f(:, 1) = h(X0, params);
newY=yNext(:,1);
xNext = X0;
X_sim(:, 1) = X0;
error=[];
tspan = [0, dt];

% Calculating the reference input for stabilizing in the reference point
uRef = mpc.computeReferenceInput();

%mhe=mhe.initialGuessPropegation();
for k=1:NT-1
    k
    if k==mhe.N_MHE+2
        mhe.Q = 5e3*mhe.Q; %This relies on having enabled dynamic update of arrival cost in MHE. 
                           %That is, G must be updated with the new Q, which is done automatically when 
                           %updating M as well in arrival cost. If this is not done, we must update G here also.
    end

    if k<=40
        [K_lqr,~,~] = dlqr(mpc.A, mpc.B, Q_LQR, R_LQR);
        %K_dlqr = [120,120,120,120,120,5,5,5,5,5;120,120,120,120,120,5,5,5,5,5;120,120,120,120,120,5,5,5,5,5;120,120,120,120,120,5,5,5,5,5];
            %Kp = K_dlqr(:,1:nStates/2);
            %Kd = K_dlqr(:,nStates/2+1:nStates);
            %U_LQR = Kp*(xEst(1:nStates/2)) + Kd*(xEst(nStates/2+1:nStates));
        U_LQR = -K_lqr*X_sim(:,k);

        U=U_LQR;

        fade_period=30:40;
        if ismember(k,fade_period)
            [~, Uopt]=mpc.runMPC(xEst);
            Udiff=U_LQR - Uopt;
            U = (1-gamma_f(k,fade_period))*U_LQR + gamma_f(k,fade_period)*Uopt;
        end
        
        [T, X] = ode45(@(t, x) f(x, U, params), tspan, X_sim(:,k));
        X_sim(:, k+1) = X(end, :)';

        U_sim(:,k) = U;
        %X_sim(:,k+1) = mpc.A*X_sim(:,k) + mpc.B*U_sim(:,k);
        newU=U_sim(:,k);

    else
        [~, Uopt]=mpc.runMPC(xEst);

        if k>=55 && k<60
            Uopt=Uopt;%+[20;20;20;20];
        end
        U_sim(:,k) = Uopt; %+ uRef;
        %X_sim(:,k+1) = mpc.A*X_sim(:,k) + mpc.B*U_sim(:,k);
        [T, X] = ode45(@(t, x) f(x, Uopt, params), tspan, X_sim(:,k));
        X_sim(:, k+1) = X(end, :)';
        newU=U_sim(:,k);
        
    end
    
    
    noise=noise_std*randn([nMeasurements,1]);
    yNext(:,k+1) = C*X_sim(:,k+1)-C*xlp+noise;
    % yNext(:, k+1) = h(X_sim(:, k+1), params);
    yNext_f(:,k+1)=alpha*yNext(:,k+1) + (1-alpha)*yNext_f(:,k);
    newY=yNext(:,k+1);
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k+1)=xEst;

end


%% Plot

figure(1)
clf
subplot(3,1,1)
plot(X_sim(1,:), "-r"); hold on; grid on;
plot(MHE_est(1,:), "-b")
ylabel("x")
title("Position")
legend(["sim","est"])

subplot(3,1,2)
plot(X_sim(2,:), "-r"); hold on; grid on;
plot(MHE_est(2,:), "-b")
ylabel("y")
legend(["sim","est"])

subplot(3,1,3)
plot(X_sim(3,:), "-r"); hold on; grid on;
plot(MHE_est(3,:)+zeq, "-b")
ylabel("z")
xlabel("iterations")
legend(["sim","est"])


figure(2)
clf
subplot(2,1,1)
plot(X_sim(4,:), "-r"); hold on; grid on;
plot(MHE_est(4,:), "-b")
ylabel("phi")
title("Angles")
legend(["sim","est"])

subplot(2,1,2)
plot(X_sim(5,:), "-r"); hold on; grid on;
plot(MHE_est(5,:), "-b")
ylabel("theta")
xlabel("iterations")
legend(["sim","est"])


figure(3)
clf
subplot(3,1,1)
plot(X_sim(6,:), "-r"); hold on; grid on;
plot(MHE_est(6,:), "-b")
ylabel("xdot")
title("Velocity")
legend(["sim","est"])

subplot(3,1,2)
plot(X_sim(7,:), "-r"); hold on; grid on;
plot(MHE_est(7,:), "-b")
ylabel("ydot")
legend(["sim","est"])

subplot(3,1,3)
plot(X_sim(8,:), "-r"); hold on; grid on;
plot(MHE_est(8,:), "-b")
ylabel("zdot")
xlabel("iterations")
legend(["sim","est"])


figure(4)
clf
subplot(2,1,1)
plot(X_sim(9,:), "-r"); hold on; grid on;
plot(MHE_est(9,:), "-b")
ylabel("phidot")
title("Angular velocity")
legend(["sim","est"])

subplot(2,1,2)
plot(X_sim(10,:), "-r"); hold on; grid on;
plot(MHE_est(10,:), "-b")
ylabel("thetadot")
xlabel("iterations")
legend(["sim","est"])


%% plot meas
mhe_meas=C*MHE_est
figure(5)
clf
subplot(3,1,1)
var=1;
plot(yNext_f(var,:)); hold on
plot(mhe_meas(var,:))
legend(["meas","est"])

subplot(3,1,2)
plot(yNext_f(var+1,:)); hold on
plot(mhe_meas(var+1,:))
legend(["meas","est"])

subplot(3,1,3)
plot(yNext_f(var+2,:)); hold on
plot(mhe_meas(var+2,:))
legend(["meas","est"])













function gamma = gamma_f(k,fade_period)
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end
