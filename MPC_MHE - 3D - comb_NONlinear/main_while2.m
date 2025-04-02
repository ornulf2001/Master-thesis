clc,clear
%addpath(genpath("3D model"))
addpath(genpath('3D model reduced order'))
addpath(genpath('../qpOASES/interfaces/matlab'))


% Define system parameters
params = parameters;

% Find equilibrium point
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1

zeq =  fzero(fz,0.1);

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

N_MHE=15;
N_MPC=5;
dt=0.003;


% States: | x y z phi theta xdot ydot zdot phidot thetadot |

alpha_NEES=0.7;
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));         
Q_MHE=10e6*diag([100,100,100,100,100,100,100,100,100,100]); 
    %Start out with low Q to trust measurements during start up, 
    %then increase Q after N_MHE+1. 
    %See below in loop
                                  

M_MHE = 1e2*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
P0 = inv(M_MHE);

Q_MPC = diag([500 500 2000 10 10 1 1 10 1 1]);
    R_MPC = diag([0.2 0.2 0.2 0.2]);

Q_LQR = diag([ ...
   1e1,1e1,1e1,1e1,1e1, ...
   1e1,1e1,1e5,1e1,1e1
   ]);
R_LQR = 1e2*eye(4);
weightScaling=1e-4;

% Bounds
run("mpc_bounds.m")

% Run
%MHE_options = qpOASES_options();
MHE_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, [], [], Q_MPC, R_MPC, nStates, nControls,MPC_options, xRef, [], []);

X_sim = zeros(nStates, NT);
U_sim = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
MHE_est(:,1)=mhe.x0; xEst = mhe.x0;
yNext=zeros(nMeasurements,NT);  
yNext(:,1)= C*(X0-xlp);
yNext_f=zeros(nMeasurements,NT);
yNext_f(:,1)=C*(X0-xlp);
NIS_traj = zeros(NT-1,1);
NEES_traj = zeros(NT-1,1);
Innovations_traj = zeros(nMeasurements,NT-1);
newY=yNext(:,1);
xNext = X0;
X_sim(:, 1) = X0;
error=[];
tspan = [0, dt];

% Calculating the reference input for stabilizing in the reference point
uRef = mpc.computeReferenceInput();

iterCounter = 1;
switchCounter = 0;
switchThreshold = 10;
NIS_current = mhe.nMeasurements;
RunningFlag = true;

dof_NIS = mhe.nMeasurements;       % degrees of freedom (number of measurements)
alpha_NIS = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NIS = chi2inv(alpha_NIS/2, dof_NIS);
upperBound_NIS = chi2inv(1 - alpha_NIS/2, dof_NIS);
useAdvancedControl = false;

[K_lqr,~,~] = dlqr(mpc.A, mpc.B, Q_LQR, R_LQR);

profile clear
profile on
while RunningFlag == true && iterCounter < (NT)
    t_start = tic;
    k=iterCounter;
    iterCounter = iterCounter + 1;
    

    if (NIS_current) >= lowerBound_NIS && (NIS_current <= upperBound_NIS)
        switchCounter = switchCounter + 1;
    else
        switchCounter = switchCounter - 3;
        if switchCounter <0
            switchCounter=0;
        end
    end

    if switchCounter > switchThreshold
        useAdvancedControl = true;
    else
        useAdvancedControl = false;
    end
    
    disp(string(k) + ", Running with advanced control: " + string(useAdvancedControl))


    if iterCounter==mhe.N_MHE+2
        mhe.Q = 5e3*mhe.Q; %This relies on having enabled dynamic update of arrival cost in MHE. 
                           %That is, G must be updated with the new Q, which is done automatically when 
                           %updating M as well in arrival cost. If this is not done, we must update G here also.
    end

    
    
    
    if useAdvancedControl
        [~, Uopt]=mpc.runMPC(xEst);
        U=Uopt;
    else
        U_LQR = -K_lqr*xEst;
        U=U_LQR;
    end
        


    [~, X] = ode45(@(t, x) f(x, U, params), tspan, X_sim(:,k));
    X_sim(:, k+1) = X(end, :)';
    U_sim(:,k) = U;
    newU=U_sim(:,k);

    noise=noise_std*randn([nMeasurements,1]);
    yNext(:,k+1) = C*X_sim(:,k+1)-C*xlp+noise;
    yNext_f(:,k+1)=alpha_NEES*yNext(:,k+1) + (1-alpha_NEES)*yNext_f(:,k);
    newY=yNext(:,k+1);
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k+1)=xEst;
    NIS_current=mhe.currentNIS;
    NIS_traj(k) = NIS_current;
    Innovations_traj(:,k)=mhe.currentInnovation;
    error=xEst - (X_sim(:,k+1)-xlp);
    NEES_traj(k) = error' / mhe.currentP * error;

    %profile viewer
    elapsed = toc(t_start)
    %RunningFlag=false;
end
profile off

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
mhe_meas=C*MHE_est;
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


%% Plot NIS chi2 

figure(6);
clf
plot(NIS_traj, 'LineWidth', 1.5); hold on;
yline(lowerBound_NIS, '--r', 'LineWidth', 1.5);
yline(upperBound_NIS, '--r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('NIS');
title(['NIS trajectory with 95% Chi-square bounds (DoF = ' num2str(mhe.nMeasurements) ')']);
grid on;
legend('NIS', 'Lower 95% bound', 'Upper 95% bound');

%% Plot NEES chi2 

dof_NEES = mhe.nStates;       % degrees of freedom (number of measurements)
alpha_NEES = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NEES = chi2inv(alpha_NEES/2, dof_NEES);
upperBound_NEES = chi2inv(1 - alpha_NEES/2, dof_NEES);

figure(7);
clf
plot(NEES_traj(35:end), 'LineWidth', 1.5); hold on;
yline(lowerBound_NEES, '--r', 'LineWidth', 1.5);
yline(upperBound_NEES, '--r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('NEES');
title(['NEES trajectory with 95% Chi-square bounds (DoF = ' num2str(mhe.nStates) ')']);
grid on;
legend('NEES', 'Lower 95% bound', 'Upper 95% bound');


%% Innovations whiteness test
%innov_var=9;
%[h,p]=lbqtest(Innovations_traj(innov_var,:));
%figure(7)
%clf
%plot(Innovations_traj(innov_var,:))
%title(['Ljung-Box test: h = '  num2str(h) ', p = ' num2str(p)])

allPassed=true;
for innov_var=1:nMeasurements
    [hj,pj]=lbqtest(Innovations_traj(innov_var,N_MHE+1:end));
    if hj~=0
        disp(['Innovations for state ' num2str(innov_var) ' did not pass the Ljung-Box test. P = ' num2str(pj)])
        allPassed=false;
        figure(8+innov_var)
        clf
        autocorr(Innovations_traj(innov_var,:), 'NumLags', 50); % z innovation
        title(num2str(innov_var))

    end
end
if allPassed==true
    disp('All states passed the Ljung-Box test! :-)')
end








function gamma = gamma_f(k,fade_period)
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end
