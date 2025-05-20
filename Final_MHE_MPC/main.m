clc,clear
%addpath(genpath("3D model"))
addpath(genpath('3D model reduced order_fixed'))
addpath(genpath('../qpOASES/interfaces/matlab'))


% Define system parameters
params = parameters1;
params_adjusted = parameters1_adjusted;

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
[Ac2, Bc2, C2] = linearizeModel(@f, @h, xlp, ulp, params_adjusted);

nStates = size(Ac, 1);
nControls = size(Bc, 2);
nMeasurements = size(C, 1);

%% Tuning
xRef = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
X0 = [0.003; 0.003; zeq+0.002; 0; 0; 0; 0; 0; 0; 0;];
MHE_x0 = zeros(nStates,1);
t = 1;

N_MHE = 15;
N_MPC = 10;
dt = 0.003;
NT = ceil(t/dt);
tvec = dt*(0:1:NT-1);
tvec2=0:dt:t+dt;

%MHE tuning
alpha = 0.9;
noise_std = 0.1 * 1e-3; %0.1 mT
R_MHE = inv(noise_std^2 * eye(nMeasurements));  
Q_MHE=1e5*diag([1,1,1,1,1,5,5,5,5,5]); 
Qscaling = 5e3;
    %Start out with low Q during start up, then increase Q after N_MHE+1. 
    %See below in loop

%Arrival cost weight initial guess (updates KF-style in loop)
%M_MHE = 1e0*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
M_MHE = 1e0*diag([5,5,5,1,1,1,1,1,1,1]);

P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling = 1e-4; %Scaling factor for better posing of QP

%MPC and LQR tuning
Q_MPC = diag([1500 1000 2000 10 10 1 1 10 1 1]);
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
mhe = MHEclass(N_MHE, Ac2, Bc2, C2, Q_MHE, R_MHE, M_MHE, weightScaling, ...
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
NIS_traj = zeros(NT, 1);
NEES_traj = zeros(NT, 1);
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
while RunningFlag == true && iterCounter <= (NT)
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
        mhe.Q = Qscaling * mhe.Q;
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
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color2='blue';
linestyle2='--';
shadecolor=[0.682, 0.847, 1];

figure(51); clf
t=tiledlayout(5,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('$\textbf{GT simulation vs. MHE estimates}$','interpreter','latex','FontSize', labelFontSize-10)

nexttile
plot(tvec2, X_sim(1,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(1,:), 'color',color2, 'lineStyle',linestyle2, 'LineWidth', lineWidth);
title('\textit{Position and orientation}','interpreter','latex','FontSize', labelFontSize-2)
set(gca,  'XTickLabel', [],'FontSize', fontSize);
ylabel('$x$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

nexttile
plot(tvec2, X_sim(2,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(2,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth', lineWidth);
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$y$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight

%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

nexttile
plot(tvec2, X_sim(3,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(3,:)+zeq, 'color',color2, 'lineStyle',linestyle2, 'LineWidth', lineWidth);
%yline(zeq, 'klinestyle2', 'LineWidth', lineWidth);
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$z$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight



nexttile
plot(tvec2,X_sim(4,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec2,MHE_est(4,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel("$\alpha$ [rad]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'location', 'northeast')
axis tight

nexttile
plot(tvec2,X_sim(5,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec2,MHE_est(5,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\beta$ [rad]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
%xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
axis tight

l = legend(t.Children(end), {'GT simulation','MHE estimates','MPC active'}, 'Orientation', 'horizontal', 'Interpreter','latex', 'FontSize', fontSize);
l.Layout.Tile = 'south'; % places legend below all plots

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 44]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 44])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch/plot_'+controllerModePrint+'_combination', '-dpdf', '-vector', '-fillpage');

%%
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color2='blue';
linestyle2='--';
shadecolor=[0.682, 0.847, 1];

figure(3); clf
t=tiledlayout(5,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('$\textbf{GT simulation vs. MHE estimates}$','interpreter','latex','FontSize', labelFontSize-10)

nexttile
plot(tvec2, X_sim(6,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(6,:), 'color',color2, 'lineStyle',linestyle2, 'LineWidth', lineWidth);
title('\textit{Translational and angular velocity}','interpreter','latex','FontSize', labelFontSize-2)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$\dot{x}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

nexttile
plot(tvec2, X_sim(7,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(7,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth', lineWidth);
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$\dot{y}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight

%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

nexttile
plot(tvec2, X_sim(8,:), 'r-', 'LineWidth', lineWidth); hold on; grid on;
plot(tvec2, MHE_est(8,:), 'color',color2, 'lineStyle',linestyle2, 'LineWidth', lineWidth);
%yline(zeq, 'klinestyle2', 'LineWidth', lineWidth);
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$\dot{z}$ [m/s]', 'Interpreter','latex', 'FontSize', labelFontSize);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2, controllModeVec,yl,shadecolor);
axis tight



nexttile
plot(tvec2,X_sim(9,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec2,MHE_est(9,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel("$\dot{\alpha}$ [rad/s]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'location', 'northeast')
axis tight

nexttile
plot(tvec2,X_sim(10,:), "r-", 'LineWidth',lineWidth); hold on; grid on;
plot(tvec2,MHE_est(10,:), 'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth)
set(gca,  'FontSize', fontSize);
ylabel("$\dot{\beta}$ [rad/s]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
%xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
%legend({'GT','MHE-based', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
axis tight

l = legend(t.Children(end), {'GT simulation','MHE estimates','MPC active'}, 'Orientation', 'horizontal', 'Interpreter','latex', 'FontSize', fontSize);
l.Layout.Tile = 'south'; % places legend below all plots

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 44]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 44])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_combination_vel', '-dpdf', '-vector', '-fillpage');


%% plot meas
%Plotting devation measurements
mhe_meas = C * MHE_est;

set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';
linestyle2='--';

figure(5);clf
t=tiledlayout(3,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\textbf{Raw measurements vs. measured estimates}','interpreter','latex','FontSize', labelFontSize)

nexttile
plot(tvec2,yNext(1,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(1,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
title('\textit{Sensor $1$}','interpreter','latex','FontSize', labelFontSize-2)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$b_x$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(2,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(2,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'XTickLabel', [],'FontSize', fontSize);
ylabel('$b_y$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(3,:),'color',color1,'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(3,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$b_z$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

l = legend(t.Children(end), {'Raw measurements','Measured estimates'}, 'Orientation', 'horizontal', 'Interpreter','latex', 'FontSize', fontSize);
l.Layout.Tile = 'south'; % places legend below all plots

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 25]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 25])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_EstMeas1', '-dpdf', '-vector', '-fillpage');

figure(6);clf
t=tiledlayout(3,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\textbf{Raw measurements vs. measured estimates}','interpreter','latex','FontSize', labelFontSize)

nexttile
plot(tvec2,yNext(4,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(4,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
title('\textit{Sensor $2$}','interpreter','latex','FontSize', labelFontSize-2)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$b_x$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(5,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(5,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'XTickLabel', [],'FontSize', fontSize);
ylabel('$b_y$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(6,:),'color',color1,'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(6,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$b_z$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

l = legend(t.Children(end), {'Raw measurements','Measured estimates'}, 'Orientation', 'horizontal', 'Interpreter','latex', 'FontSize', fontSize);
l.Layout.Tile = 'south'; % places legend below all plots

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 25]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 25])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_EstMeas2', '-dpdf', '-vector', '-fillpage');

figure(7);clf
t=tiledlayout(3,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\textbf{Raw measurements vs. measured estimates}','interpreter','latex','FontSize', labelFontSize)

nexttile
plot(tvec2,yNext(7,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(7,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
title('\textit{Sensor $3$}','interpreter','latex','FontSize', labelFontSize-2)
set(gca, 'XTickLabel', [], 'FontSize', fontSize);
ylabel('$b_x$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(8,:),'color',color1, 'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(8,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'XTickLabel', [],'FontSize', fontSize);
ylabel('$b_y$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

nexttile
plot(tvec2,yNext(9,:),'color',color1,'LineWidth',lineWidth-0.5); hold on; grid on
plot(tvec2,mhe_meas(9,:),'color',color2, 'lineStyle',linestyle2,'LineWidth',lineWidth);
set(gca,  'FontSize', fontSize);
ylabel('$b_z$ [T]', 'Interpreter','latex', 'FontSize', labelFontSize);
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

l = legend(t.Children(end), {'Raw measurements','Measured estimates'}, 'Orientation', 'horizontal', 'Interpreter','latex', 'FontSize', fontSize);
l.Layout.Tile = 'south'; % places legend below all plots

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 25]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 25])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_EstMeas3', '-dpdf', '-vector', '-fillpage');



%meas_error=vecnorm(yNext - mhe_meas);
%plot(tvec,meas_error(:,:))
%title("Vector norm of Y-C*xEst")

%% Plot NIS chi2 
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 23;
labelFontSize = 23;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';
linestyle2='--';
shadecolor=[0.682, 0.847, 1];
cap = 2 * upperBound_NIS;

figure(6);
clf

t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
titleStr = '$\textbf{NIS trajectory with 95\% }\chi^2 \textbf{ bounds (DoF = }$'+ string(mhe.nMeasurements)+')';
sgtitle(titleStr, 'Interpreter', 'latex', 'FontSize', labelFontSize);

nexttile
plot(tvec,NIS_traj, 'color',color2, 'LineWidth', lineWidth); hold on; grid on;
yline(lowerBound_NIS, '--r', 'LineWidth', lineWidth);
yline(upperBound_NIS, '--r', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
ylabel('NIS', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

ylim([-0,cap])
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
for k = 1:length(NIS_traj)
    if NIS_traj(k) > cap
        plot(tvec(k), 0.994*cap, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % upside-down triangle
    end
end

legend({'NIS', '95\% bounds','',"MPC active",'','','','','','','','NIS above axis limit'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_NIS', '-dpdf', '-vector', '-fillpage');

%% Plot NEES chi2 
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 23;
labelFontSize = 23;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';
linestyle2='--';
shadecolor=[0.682, 0.847, 1];
cap = 2 * upperBound_NIS;

dof_NEES = mhe.nStates;       % degrees of freedom (number of measurements)
alpha_NEES = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NEES = chi2inv(alpha_NEES/2, dof_NEES);
upperBound_NEES = chi2inv(1 - alpha_NEES/2, dof_NEES);

figure(6);
clf

t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
titleStr = '$\textbf{NEES trajectory with 95\% }\chi^2 \textbf{ bounds (DoF = }$'+ string(mhe.nStates)+')';
sgtitle(titleStr, 'Interpreter', 'latex', 'FontSize', labelFontSize);

nexttile
plot(tvec,NEES_traj, 'color',color2, 'LineWidth', lineWidth); hold on; grid on;
yline(lowerBound_NEES, '--r', 'LineWidth', lineWidth);
yline(upperBound_NEES, '--r', 'LineWidth', lineWidth);
set(gca,  'FontSize', fontSize);
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize);
ylabel('NEES', 'Interpreter','latex', 'FontSize', labelFontSize);
axis tight

ylim([-0,cap])
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec2,controllModeVec,yl,shadecolor)
for k = 1:length(NEES_traj)
    if NEES_traj(k) > cap
        plot(tvec(k), 0.994*cap, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % upside-down triangle
    end
end

legend({'NEES', '95\% bounds','',"MPC active",'','','','','','','NEES above axis limit'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_NEES', '-dpdf', '-vector', '-fillpage');

% dof_NEES = mhe.nStates;       % degrees of freedom (number of measurements)
% alpha_NEES = 0.05;  % 95% confidence = 1 - alpha
% lowerBound_NEES = chi2inv(alpha_NEES/2, dof_NEES);
% upperBound_NEES = chi2inv(1 - alpha_NEES/2, dof_NEES);
% 
% figure(7);
% clf
% plot(NEES_traj(:), 'color',color2, 'LineWidth', lineWidth); hold on;grid on;
% yline(lowerBound_NEES, '--r', 'LineWidth', lineWidth);
% yline(upperBound_NEES, '--r', 'LineWidth', lineWidth);
% xlabel('Time', 'Interpreter','latex', 'FontSize', labelFontSize);
% ylabel('NEES', 'Interpreter','latex', 'FontSize', labelFontSize);
% title(['NEES trajectory with 95\% Chi-square bounds (DoF = ' num2str(mhe.nStates) ')'], 'FontSize',labelFontSize+2, 'Interpreter','latex');
% ylim([-0,70])
% yl = get(gca, 'YLim');ylim(yl)
% shadeMPCregions(tvec,controllModeVec,yl)
% legend({'NEES', 'Lower 95\% bound', 'Upper 95\% bound', 'MPC active'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast');
% print(gcf, 'Figures/estimation_error/plot_'+controllerModePrint+'_NEES', '-dsvg');


%% Innovations whiteness test

set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 23;
labelFontSize = 23;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';

allPassed=true;
for innov_var=1:nMeasurements
    [hj,pj]=lbqtest(Innovations_traj(innov_var,ceil(0.2/dt):end));
    if hj~=0
        disp(['Innovations for measurement ' num2str(innov_var) ' did not pass the Ljung-Box test. P = ' num2str(pj)])
        allPassed=false;

        figure(8+innov_var)
        clf
        t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
        titleStr = '$\textbf{ACF of innovations for measurement '+ string(innov_var)+' (Ljung-Box failed)}$';
        sgtitle(titleStr, 'Interpreter', 'latex', 'FontSize', labelFontSize);

        autocorr(Innovations_traj(innov_var,:), 'NumLags', 50); % z innovation
        set(gca,  'FontSize', fontSize);
        xlabel('Lag', 'Interpreter','latex', 'FontSize', labelFontSize);
        ylabel('Autocorrelation', 'Interpreter','latex', 'FontSize', labelFontSize);
        %title(['Measurement ' num2str(innov_var)], 'FontName', 'Times New Roman', 'FontSize', labelFontSize+2);
        title('')
        set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
        set(gcf, 'PaperUnits', 'centimeters')
        set(gcf, 'PaperSize', [41 20])
        set(gcf, 'PaperPositionMode', 'manual')
        print(gcf, 'Figures/model_mismatch//plot_'+controllerModePrint+'_ACF'+string(innov_var),'-dpdf', '-vector', '-fillpage');

    end
end
if allPassed==true
    disp('All innovations passed the Ljung-Box test! :-)')
end


%%
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 23;
labelFontSize = 23;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
shadecolor=[0.682, 0.847, 1];

color2='blue';

t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
titleStr = '$\textbf{Control inputs for GT simulation}$';
sgtitle(titleStr, 'Interpreter', 'latex', 'FontSize', labelFontSize);

nexttile
figure(100);clf
plot(tvec,U_sim(:,:), 'LineWidth', lineWidth); hold on
set(gca,  'FontSize', fontSize);
ylabel("Solenoid current [A]", 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlabel("Iterations", 'Interpreter','latex', 'FontSize', labelFontSize)
%title("Control input", 'FontSize',labelFontSize+2)
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl,shadecolor)
legend({"$u_{x,p}$","$u_{y,p}$","$u_{x,n}$","$u_{y,n}$","MPC active"}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast')

%Zoom-in box (Manually located)
axInset = axes('Position', [0.3, 0.25, 0.1, 0.5]);  % [x y width height]
box on;
plot(tvec, U_sim(:,:), 'LineWidth', lineWidth);
yl = get(gca, 'YLim');ylim(yl)
shadeMPCregions(tvec,controllModeVec,yl,shadecolor)
xlim([0.1, 0.2]);   % zoomed x-range
ylim([-6, 1]);   % zoomed y-range
set(axInset, 'FontSize', 10);

% Draw connecting lines (normalized coordinates)
annotation('line', [0.3 0.225], [0.50 0.65], 'LineStyle', '--');
annotation('line', [0.3  0.225], [0.75 0.65], 'LineStyle', '--');

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/estimation_error/plot_'+controllerModePrint+'_control', '-dpdf', '-vector', '-fillpage');


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