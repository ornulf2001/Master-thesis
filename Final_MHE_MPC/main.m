clc,clear
%addpath(genpath("3D model"))
addpath(genpath('3D model reduced order'))
addpath(genpath('../qpOASES/interfaces/matlab'))


% Define system parameters
params = parameters;

% Find equilibrium point
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]', [0,0,0,0]', params), 8);  
zeq =  fzero(fz,0.1);

xeq = [0, 0, zeq, zeros(1,7)]';
ueq = [0,0,0,0]';

% Linearize model
xlp = xeq;   
ulp = ueq;

% States: [ x y z phi theta xdot ydot zdot phidot thetadot ]
[Ac, Bc, C] = linearizeModel(@f, @h, xlp, ulp, params);

nStates = size(Ac, 1);
nControls = size(Bc, 2);
nMeasurements = size(C, 1);

%% Tuning
xRef = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
X0 = [0.003; 0.003; zeq+0.002; 0; 0; 0; 0; 0; 0; 0;];
MHE_x0 = zeros(nStates,1);
t = 2;

N_MHE = 10;
N_MPC = 10;
dt = 0.003;
NT = 350;
tvec = 0:1:NT-1;

%MHE tuning
alpha = 0.9;
noise_std = 0.1 * 1e-3; %0.1 mT
R_MHE = inv(noise_std^2 * eye(nMeasurements));  
Q_MHE=10e3*diag([100,100,100,100,100,500,500,500,500,500]); 
    %Start out with low Q during start up, then increase Q after N_MHE+1. 
    %See below in loop

%Arrival cost weight initial guess (updates KF-style in loop)
M_MHE = 1e0*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling = 1e-4; %Scaling factor for better posing of QP

%MPC and LQR tuning
Q_MPC = diag([500 500 2000 10 10 1 1 10 1 1]);
R_MPC = diag([0.2, 0.2, 0.2, 0.2]);

Q_LQR = diag([ ...
   1e1,1e1,1e1,1e1,1e1, ...
   1e1,1e1,1e5,1e1,1e1
   ]);
R_LQR = 1e2 * eye(4);

% Bounds
run("mpc_bounds.m") %currently inf all over

%% Run

%MHE_options = qpOASES_options();
MHE_options = optimset('Display','off', 'Diagnostics','off', ...
        'Algorithm', 'active-set');
mhe = MHEclass(N_MHE, Ac, Bc, C, Q_MHE, R_MHE, M_MHE, weightScaling, ...
        MHE_x0, xlp, P0, dt, MHE_options);

MPC_options = optimset('Display', 'off', 'Diagnostics', 'off', ...
        'Algorithm', 'active-set');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, [], [], Q_MPC, R_MPC, ...
        nStates, nControls, MPC_options, xRef, [], []);

%Init
X_sim = zeros(nStates, NT);
U_sim = zeros(nControls, NT-1);
vsol = zeros(nMeasurements, NT);
wsol = zeros(nMeasurements, NT-1);
MHE_est = zeros(nStates, NT);
MHE_est(:,1) = mhe.x0; 
xEst = mhe.x0;
yNext = zeros(nMeasurements, NT);  
yNext(:,1) = C * (X0-xlp);
yNext_f = zeros(nMeasurements, NT);
yNext_f(:,1) = C * (X0-xlp);
NIS_traj = zeros(NT-1, 1);
NEES_traj = zeros(NT-1, 1);
Innovations_traj = zeros(nMeasurements, NT-1);
newY = yNext(:,1);
xNext = X0;
X_sim(:, 1) = X0;
tspan = [0, dt];
controllModeVec = zeros(1,NT);
toterror=0;

% Calculating the reference input for stabilizing in the reference point
uRef = mpc.computeReferenceInput(); 

iterCounter = 1;
switchCounter = 0;
switchThreshold = 10;
NIS_current = mhe.nMeasurements;
RunningFlag = true;

dof_NIS = mhe.nMeasurements; % degrees of freedom (number of measurements)
alpha_NIS = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NIS = chi2inv(alpha_NIS / 2, dof_NIS);
upperBound_NIS = chi2inv(1 - alpha_NIS / 2, dof_NIS);
useAdvancedControl = false;

[K_lqr,~,~] = dlqr(mpc.A, mpc.B, Q_LQR, R_LQR);

profile clear
profile on

% Current switching logic: 
%   If NIS is inside bounds: increase switchCounter, if Nis exits bounds: reduce switchCounter. 
%   Whenever switchCounter > switchThreshold: use advanced control, else use LQR.
%   This hopefully leads to a smoother transition between control types.
while RunningFlag == true && iterCounter < (NT)
    t_start = tic;
    k = iterCounter;
    iterCounter = iterCounter + 1;
    

    if (NIS_current) >= lowerBound_NIS && (NIS_current <= upperBound_NIS)
        switchCounter = switchCounter + 1;
        if switchCounter > 2 * switchThreshold
            switchCounter = 2 * switchThreshold;
        end
    else
        switchCounter = switchCounter - 3;
        if switchCounter < 0
            switchCounter = 0;
        end
    end

    if switchCounter > switchThreshold
        useAdvancedControl = true;
    else
        useAdvancedControl = false;
    end
    
    disp(string(k) + ", Running with advanced control: " + ...
            string(useAdvancedControl))


    if iterCounter == mhe.N + 2
        mhe.Q = 5e3 * mhe.Q;
        mhe.G(mhe.nStates * (mhe.N+1) + 1 : mhe.nStates * (mhe.N+1) + ...
                mhe.nStates * mhe.N, mhe.nStates * (mhe.N+1) + 1 : ...
                mhe.nStates * (mhe.N+1) + mhe.nStates * mhe.N ) = ...
                kron(eye(mhe.N), mhe.weightScaling * mhe.Q);
        % Here we increase Q after N_MHE+1 iterations when the MHE has calibrated. 
        % This seems to improve performance?
    end

    %use controllerMode = X_sim(:,k) or xEst for running MHE on true state or MHE estimates
    controllerMode = X_sim(:,k); if controllerMode ==xEst; controllerModePrint="xEst"; elseif controllerMode == X_sim(:,k); controllerModePrint="Xsim";end
    if useAdvancedControl
        [~, Uopt] = mpc.runMPC(controllerMode);
        U = Uopt;
        controllModeVec(k) = 1;
    else
        U_LQR = -K_lqr * controllerMode;
        U = U_LQR;
    end
        


    %[~, X] = ode15s(@(t, x) f(x, U, params), tspan, X_sim(:,k));
    %X_sim(:, k+1) = X(end, :)'; 
    X_sim(:, k+1) = RK4Step(@f, X_sim(:,k), U, dt, params);
    U_sim(:, k) = U; 
    newU = U_sim(:, k); %For MHE input

    noise = noise_std * randn([nMeasurements, 1]);
    yNext(:, k+1) = C * X_sim(:, k+1) - C * xlp + noise; 
    yNext_f(:, k+1) = alpha * yNext(:, k+1) + (1-alpha) * yNext_f(:, k); %EMA
    newY = yNext_f(:, k+1); %For MHE input
    mhe = mhe.runMHE(newY, newU); %Run MHE with newY and newU
    xEst = mhe.xCurrent; % Extract xhat
    MHE_est(:, k+1) = xEst;
    vsol(:, k+1) = mhe.vCurrent;
    wsol(:, k+1) = mhe.vCurrent;


    NIS_current = mhe.currentNIS;
    NIS_traj(k) = NIS_current;
    Innovations_traj(:,k) = mhe.currentInnovation;
    error = xEst - (X_sim(:, k+1) - xlp);
    toterror = toterror + vecnorm(error,2);
    NEES_traj(k) = error' / mhe.currentP * error;

    %profile viewer
    elapsed = toc(t_start)
end
%%
error_traj = (X_sim-xlp)-MHE_est;
RMSE = sqrt(mean(error_traj.^2,2));
RMSE2 = sqrt(mean(error(:).^2));
profile off
save("data/Y_noisy_sim","yNext_f")
save("data/U_list_sim","U_sim")

%% Plot
% Shared settings
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 14;
labelFontSize = 18;
lineWidth = 1.5;


figure(1); clf

% --- Subplot 1: x ---
subplot(3,1,1)
plot(tvec, X_sim(1,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(1,:), 'b--', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$x$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
title('Position', 'FontSize', labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec, controllModeVec,yl);
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Subplot 2: y ---
subplot(3,1,2)
plot(tvec, X_sim(2,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(2,:), 'b--', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$y$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
shadeMPCregions(tvec, controllModeVec,yl);
yl = get(gca, 'YLim');ylim(yl)
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Subplot 3: z ---
subplot(3,1,3)
plot(tvec, X_sim(3,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(3,:) + zeq, 'b--', 'LineWidth', lineWidth);
yline(zeq, 'k-.', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$z$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
xlabel('Iterations', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec, controllModeVec,yl);
legend({'Simulation','Estimates','$\texttt{zeq}$', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Export ---
print(gcf, 'Figures/plot_'+controllerModePrint+'_Position', '-dsvg');


%%
figure(2)
clf
subplot(2,1,1)
plot(tvec,X_sim(4,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec,MHE_est(4,:), "b--",'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\alpha$ [rad]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
title("Angle", 'FontSize',labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({"Simulation","Estimates", 'MPC active'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'location', 'northeast')

subplot(2,1,2)
plot(tvec,X_sim(5,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec,MHE_est(5,:), "b--",'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\beta$ [rad]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

print(gcf, 'Figures/plot_'+controllerModePrint+'_Angle', '-dsvg');

%%

figure(3); clf

% --- Subplot 1: x ---
subplot(3,1,1)
plot(tvec, X_sim(6,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(6,:), 'b--', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$\dot{x}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
title('Velocity', 'FontSize', labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec, controllModeVec,yl);
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Subplot 2: y ---
subplot(3,1,2)
plot(tvec, X_sim(7,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(7,:), 'b--', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$\dot{y}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
shadeMPCregions(tvec, controllModeVec,yl);
yl = get(gca, 'YLim');ylim(yl)
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Subplot 3: z ---
subplot(3,1,3)
plot(tvec, X_sim(8,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec, MHE_est(8,:), 'b--', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$\dot{z}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
xlabel('Iterations', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec, controllModeVec,yl);
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

% --- Export ---
print(gcf, 'Figures/plot_'+controllerModePrint+'_Velocity', '-dsvg');




figure(4)
clf
subplot(2,1,1)
plot(tvec,X_sim(9,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec,MHE_est(9,:), "b--",'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\dot{\alpha}$ [rad/s]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
title("Angular velocity", 'FontSize',labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({"Simulation","Estimates", 'MPC active'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'location', 'northeast')

subplot(2,1,2)
plot(tvec,X_sim(10,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec,MHE_est(10,:), "b--",'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\dot{\beta}$ [rad/s]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({'Simulation','Estimates', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
print(gcf, 'Figures/plot_'+controllerModePrint+'_AngularVelocity', '-dsvg');


%% plot meas
%Plotting devation measurements
mhe_meas=C*MHE_est;
meas_error=vecnorm(yNext_f-mhe_meas,2);
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

figure(101);clf
%plot(tvec,yNext_f(1:9,:),"r-"); hold on
%plot(tvec,mhe_meas(1:9,:),"b--")
plot(tvec,meas_error(:,:))
title("Vector norm of Y-C*xEst")
print(gcf, 'Figures/plot_'+controllerModePrint+'_EstMeas', '-dsvg');


%% Plot NIS chi2 

figure(6);
clf
plot(NIS_traj, 'LineWidth', lineWidth); hold on;
yline(lowerBound_NIS, '--r', 'LineWidth', lineWidth);
yline(upperBound_NIS, '--r', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
xlabel('Time', 'Interpreter','latex', 'FontSize', labelFontSize);
ylabel('NIS', 'Interpreter','latex', 'FontSize', labelFontSize);
title(['NIS trajectory with 95\% Chi-square bounds (DoF = ' num2str(mhe.nMeasurements) ')'], ...
    'FontSize', labelFontSize+2, 'Interpreter', 'latex');
grid on;
ylim([-0,40])
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({'NIS', 'Lower 95\% bound', 'Upper 95\% bound', "MPC active"}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
print(gcf, 'Figures/plot_'+controllerModePrint+'_NIS', '-dsvg');

%% Plot NEES chi2 

dof_NEES = mhe.nStates;       % degrees of freedom (number of measurements)
alpha_NEES = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NEES = chi2inv(alpha_NEES/2, dof_NEES);
upperBound_NEES = chi2inv(1 - alpha_NEES/2, dof_NEES);

figure(7);
clf
plot(NEES_traj(:), 'LineWidth', lineWidth); hold on;grid on;
yline(lowerBound_NEES, '--r', 'LineWidth', lineWidth);
yline(upperBound_NEES, '--r', 'LineWidth', lineWidth);
xlabel('Time', 'Interpreter','latex', 'FontSize', labelFontSize);
ylabel('NEES', 'Interpreter','latex', 'FontSize', labelFontSize);
title(['NEES trajectory with 95\% Chi-square bounds (DoF = ' num2str(mhe.nStates) ')'], 'FontSize',labelFontSize+2, 'Interpreter','latex');
ylim([-0,70])
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({'NEES', 'Lower 95\% bound', 'Upper 95\% bound', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
print(gcf, 'Figures/plot_'+controllerModePrint+'_NEES', '-dsvg');


%% Innovations whiteness test

allPassed=true;
for innov_var=1:nMeasurements
    [hj,pj]=lbqtest(Innovations_traj(innov_var,60:end));
    if hj~=0
        disp(['Innovations for measurement ' num2str(innov_var) ' did not pass the Ljung-Box test. P = ' num2str(pj)])
        allPassed=false;
        figure(8+innov_var)
        clf
        autocorr(Innovations_traj(innov_var,:), 'NumLags', 50); % z innovation
        ax = gca;
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize);
        xlabel('Lag', 'FontName', 'Times New Roman', 'FontSize', labelFontSize);
        ylabel('Autocorrelation', 'FontName', 'Times New Roman', 'FontSize', labelFontSize);
        title(['Measurement ' num2str(innov_var)], 'FontName', 'Times New Roman', 'FontSize', labelFontSize+2);
        print(gcf, 'Figures/plot_'+controllerModePrint+'_ACF'+string(innov_var), '-dsvg');

    end
end
if allPassed==true
    disp('All innovations passed the Ljung-Box test! :-)')
end


%%
figure(100);clf
plot(tvec(1:end-1),U_sim(:,:), 'LineWidth', lineWidth); hold on
set(gca,  'FontSize', fontSize);
ylabel("Solenoid current [A]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
title("Control input", 'FontSize',labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
legend({"$u_{x,p}$","$u_{y,p}$","$u_{x,n}$","$u_{y,n}$","MPC active"}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast')

%Zoom-in box (Manually located)
axInset = axes('Position', [0.3, 0.55, 0.25, 0.25]);  % [x y width height]
box on;
plot(tvec(1:end-1), U_sim(:,:), 'LineWidth', lineWidth);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl)
xlim([25, 55]);   % zoomed x-range
ylim([-0.8, 1.2]);   % zoomed y-range
set(axInset, 'FontSize', 10);

% Draw connecting lines (normalized coordinates)
annotation('line', [0.3 0.21], [0.55 0.42], 'LineStyle', '--');
annotation('line', [0.3+0.25  0.21], [0.55 0.42], 'LineStyle', '--');



%%

function gamma = gamma_f(k,fade_period) %Not in use, old fading factor [0,1] to switch from LQR to MPC
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end

function x_next = RK4Step(f,x,U,dt,params)
    k1 = f(x,U,params);
    k2 = f(x+0.5*dt*k1,U,params);
    k3 = f(x+0.5*dt*k2,U,params);
    k4 = f(x + dt*k3, U, params);

    x_next = x+(dt/6)*(k1+2*k2+2*k3+k4);
end