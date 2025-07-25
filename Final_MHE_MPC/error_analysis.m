clc, clear
addpath(genpath('3D model reduced order_fixed'))

N = 4:2:30;           % Range of N
dt = 0.003;%0.001:0.001:0.003;          % Range of dt (can be vector later)


[NGrid, dtGrid] = meshgrid(N, dt);
[nDt, nN] = size(NGrid);
nSim = 5;            % Number of Monte Carlo simulations

    NT = 500;
    % Define system parameters
    params = parameters;
    
    % Find equilibrium point
    index = @(A,i) A(i);
    fz = @(z) index(f([0,0,z,zeros(1,7)]', [0,0,0,0]', params), 8);  
    zeq =  fzero(fz,0.1);
    xeq = [0, 0, zeq, zeros(1,7)]'; xlp=xeq;
    X0 = [0.003; 0.003; zeq+0.002; 0; 0; 0; 0; 0; 0; 0;];
    MHE_x0 = X0-xlp;
    MHE_x0 = MHE_x0+[0.001;-0.001;0.001;zeros(7,1)];

RMSE = zeros(nDt, nN, nSim);
sims = zeros(10, NT, numel(NGrid),nSim);
ests = zeros(10, NT, numel(NGrid),nSim);
runtimes = zeros(nDt,nN,nSim);
%%
for k = 1:nSim
    for i = 1:numel(NGrid)
        [row, col] = ind2sub(size(NGrid), i);
        currentN = NGrid(row, col);
        currentDt = dtGrid(row, col);

        disp("Simulation series " + k)
        disp("Running simulation with N = " + ...
                num2str(currentN) + " and dt = " + num2str(currentDt))

        tic
        [sim, est] = runSystem(currentN, currentDt, NT, X0, MHE_x0,zeq);
        sims(:,:,i,k) = sim;
        ests(:,:,i,k) = est;

        %error = sim(:,NT/2:NT) - est(:,NT/2:NT);
        error = sim - (est+xlp);
        RMSE(row, col, k) = sqrt(mean(error(:).^2));
        runtimes(row,col,k)=toc;
    end
end

disp("Finished all simulations")
%%
% Average across Monte Carlo runs
RMSE_mean = mean(RMSE, 3)';
RMSE_std = std(squeeze(RMSE), 0 ,2);
%%
    tvec=dt*(1:NT+dt);

%% load previous data

% clc,clear
% datafolder = "data_20250510_140018";
% load('data/'+datafolder+'/sims.mat')
% load('data/'+datafolder+'/runtimes.mat')
% load('data/'+datafolder+'/RMSEtotal.mat')     
% load('data/'+datafolder+'/RMSEmean.mat')
% load('data/'+datafolder+'/ests.mat')

%%
folder="data/data_" + datestr(datetime("now"), "yyyymmdd_HHMMSS");
if ~exist(folder, 'dir')
    mkdir(folder);
end
save(folder + "/RMSEtotal","RMSE")
save(folder + "/sims","sims")
save(folder + "/ests","ests")
save(folder + "/RMSEmean","RMSE_mean")
save(folder + '/runtimes','runtimes')
save(folder + "/Ns","N")
save(folder+ "/dts", "dt")
save(folder+ "/MCs","nSim")
%% plot

set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color2='blue';
linestyle2='--';

figure(1)
clf
t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('$\textbf{RMSE vs. Horizon length $\texttt{N}$.}$','interpreter','latex','FontSize', labelFontSize+5)

nexttile
% try %Generally fails of length(dt)=1 -> RMSE is vector not grid
%     surf(NGrid, dtGrid, RMSE_mean); grid on; hold on
%    % surf(NGrid, dtGrid, mean(runtimes,3))
%     colorbar
%     xlabel('Horizon Length N');
%     ylabel('Time Step dt');
%     zlabel('RMSE');
%     title('Surface plot of RMSE vs. N and dt');
% 
% catch
    yyaxis left
    hold on;
    x = [N, fliplr(N)];
    y = [RMSE_mean + RMSE_std; flipud(RMSE_mean - RMSE_std)]';
    fill(x, y, [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);


    plot(N,RMSE_mean,'LineWidth', lineWidth); grid on; 
    set(gca,'FontSize', fontSize);
    ylabel('RMSE', 'Interpreter','latex', 'FontSize', labelFontSize);


    yyaxis right
    plot(N,mean(squeeze(runtimes),2),'LineWidth', lineWidth)
    ylabel("Runtimes [s]", 'Interpreter','latex', 'FontSize', labelFontSize)
    xlabel('Horizon Length $\texttt{N}$', 'Interpreter','latex', 'FontSize', labelFontSize);
    axis tight

    legend({"RMSE $\pm 1\sigma$", 'RMSE mean',"Runtimes"}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast' )
    
set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/RMSE/mid_R_osc', '-dpdf', '-vector', '-fillpage');

savefig(folder + '/RMSE_vs_N_mid_R_osc.fig')

%%
% Only use dt == dt(end) for each N (assuming dt is scalar, otherwise filter for it)
nStates = size(sims,1);
RMSE_states = zeros(length(N), nStates);
for nIdx = 1:length(N)
    % Find the correct index for this N and dt
    idx = find(NGrid(:) == N(nIdx) & dtGrid(:) == dt(end));
    idx = idx(1); % Safety, should only be one
    err_state_allMC = zeros(nStates, NT, nSim);
    for k = 1:nSim
        sim_this = sims(:,:,idx,k);  % [nStates, NT]
        est_this = ests(:,:,idx,k);  % [nStates, NT]
        err = sim_this - (est_this + xlp); % [nStates, NT]
        err_state_allMC(:,:,k) = err;
    end
    % RMSE for each state, averaged over time and MC
    RMSE_states(nIdx,:) = sqrt(mean(mean(err_state_allMC.^2,3),2))';
end

RMSE_states_norm = RMSE_states ./ max(RMSE_states, [], 1);

% (Optional) Normalize state names for plotting
stateNames = {'$x$', '$y$', '$z$', '$\alpha$', '$\beta$', '$\dot{x}$', '$\dot{y}$', '$\dot{z}$', '$\dot{\alpha}$', '$\dot{\beta}$'};
baseColors = [
    0.2, 0.4, 0.8;    % x/xdot (blue)
    0.85, 0.33, 0.1;  % y/ydot (orange)
    0.56, 0.345, 1;   % z/zdot (purple)
    0.47, 0.67, 0.19; % phi/phidot (green)
    0.94, 0.83, 0.07  % theta/thetadot (yellow)
];
idx_min = 1;
idx_mid = round(length(N)/2);
idx_max = length(N);
N_selected = [N(idx_min), N(idx_mid), N(idx_max)];
RMSE_states_norm_selected = RMSE_states_norm([idx_min, idx_mid, idx_max], :);

barColors = baseColors;

% Velocity colors: lighter version (mix with white)
velColors = 0.5*baseColors + 0.5;   % lighter shades

% Combine: first 5 pos, last 5 vel
allColors = [barColors; velColors];

% Plot grouped bar
figure(10); clf;

t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('$\textbf{Trends of normalized per-state RMSE for selected $\texttt{N}$}$','interpreter','latex','FontSize', labelFontSize+5)

nexttile

set(gcf, 'Units', 'centimeters', 'Position', [0 0 35 16])
b = bar(N_selected, RMSE_states_norm_selected, 'grouped');
hold on

for i = 1:10
    if i <= 5
        b(i).FaceColor = barColors(i,:);
        b(i).FaceAlpha = 1.0;
    else
        b(i).FaceColor = velColors(i-5,:);
        b(i).FaceAlpha = 0.7;  % Can use <1 for extra effect
    end
end

% Legend, axes, etc
xlabel('Horizon Length $\texttt{N}$', 'Interpreter', 'latex', 'FontSize', labelFontSize)
ylabel('Normalized RMSE', 'Interpreter', 'latex', 'FontSize', labelFontSize)
legend(stateNames, 'Interpreter', 'latex', 'FontSize', fontSize, 'Location', 'northeastoutside')
set(gca, 'FontSize', fontSize)
%title('Normalized RMSE per state (selected $N$)', 'Interpreter', 'latex', 'FontSize', labelFontSize+5)
grid on

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/RMSE/barplot_rmse_per_state_vs_N_normalized_selected_grouped', '-dpdf', '-vector', '-fillpage');

savefig(folder + '/barplot_rmse_per_state_vs_N_normalized_selected_grouped.fig')

%%
z_traj_all = squeeze(ests(3, :, :, :)); % [NT, numel(NGrid), nSim]
z_mean_by_N = zeros(length(N), NT);

for nIdx = 1:length(N)
    currentN = N(nIdx);

    % Use only the first dt (for each N)
    idx = find(NGrid(:) == currentN & dtGrid(:) == dt(end));

    % If for some reason multiple indices, just pick the first (should only be one)
    idx = idx(1);

    % Extract z-trajectories for all MC runs at this N and dt
    z_data = squeeze(ests(3, :, nIdx, :));   % [NT, nSim]

    % Average over MC runs
    z_mean_by_N(nIdx, :) = mean(z_data, 2); % [NT, 1]
end

sim_idx_N = find(dtGrid(:) == dt(end));  % These are indices for all N at dt(1)

%k = 1;  % Or loop/average over k if you want MC mean too
sims_meanMC=mean(sims,4);
z_sim_trajectories = zeros(length(sim_idx_N), NT);

for i = 1:length(sim_idx_N)
    idx = sim_idx_N(i);
    z_sim_trajectories(i, :) = sims_meanMC(3, :, idx);  % [1, NT]
end
sim_mean = mean(sims, 4);  % [1, NT]


idx_min = 1;
idx_mid = round(length(N)/2);
idx_max = length(N);

%%
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';
color3=[0.329, 0.329, 1];
color4=[0.56, 0.345, 1];
color5=[0.56, 1, 0.41];
linestyle2='--';

figure(2); clf
t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('$\textbf{Comparison of $z$ trajectories for selected $N$}$','interpreter','latex','FontSize', labelFontSize+5)

plot(tvec,z_mean_by_N(idx_min,:)+zeq, 'color','m', 'LineWidth', lineWidth);hold on;grid on
plot(tvec,z_mean_by_N(idx_mid,:)+zeq, 'color',color2, 'LineWidth', lineWidth)
plot(tvec,z_mean_by_N(idx_max,:)+zeq, 'color',color5, 'LineWidth', lineWidth)
yline(zeq, 'k--')
plot(tvec,z_sim_trajectories(idx_min,:), 'r-', 'LineWidth', lineWidth-0.5)
xlim([-20*dt,NT*dt])
legend({ ['N=' num2str(N(idx_min))], ...
               ['N=' num2str(N(idx_mid))], ...
               ['N=' num2str(N(idx_max))], 'zeq','sim'}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast')
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize); 
ylabel('$z$ [m]', 'Interpreter','latex', 'FontSize', labelFontSize)
set(gca,'FontSize', fontSize);

%Zoom-in box (Manually located)
axInset = axes('Position', [0.58, 0.22, 0.3, 0.15]);  % [x y width height]
box on;
plot(tvec,z_mean_by_N(idx_min,:)+zeq, 'color','m', 'LineWidth', lineWidth);hold on;grid on
plot(tvec,z_mean_by_N(idx_mid,:)+zeq, 'color',color2, 'LineWidth', lineWidth)
plot(tvec,z_mean_by_N(idx_max,:)+zeq, 'color',color5, 'LineWidth', lineWidth)
plot(tvec,z_sim_trajectories(idx_min,:), 'r-', 'LineWidth', lineWidth)


xlim([1.3, 1.47]);   % zoomed x-range
ylim([0.0317, 0.0319]);   % zoomed y-range
set(axInset, 'FontSize', fontSize);

% Draw connecting lines (normalized coordinates)
annotation('line', [0.58 0.85], [0.37 0.48], 'LineStyle', '--');
annotation('line', [0.88  0.85], [0.37 0.48], 'LineStyle', '--');

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/RMSE/mid_R_osc_ztraj', '-dpdf', '-vector', '-fillpage');

savefig(folder + '/ztraj_mid_R_osc_selectedN.fig')
%%




set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fontSize = 22;
labelFontSize = 22;
lineWidth = 3;
color1=[0.85, 0.33, 0.1];
color2='blue';
color3=[0.329, 0.329, 1];
color4=[0.56, 0.345, 1];
color5=[0.56, 1, 0.41];
linestyle2='--';

figure(3); clf
t=tiledlayout(1,1,'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle('\textbf{Estimation Error Trajectories}','interpreter','latex','FontSize', labelFontSize+5)

plot(tvec,z_mean_by_N(idx_min,:)+zeq-z_sim_trajectories(idx_min,:), 'color','m', 'LineWidth', lineWidth);hold on;grid on
plot(tvec,z_mean_by_N(idx_mid,:)+zeq-z_sim_trajectories(idx_mid,:), 'color',color2, 'LineWidth', lineWidth)
plot(tvec,z_mean_by_N(idx_max,:)+zeq-z_sim_trajectories(idx_max,:), 'color',color5, 'LineWidth', lineWidth)
xlim([-20*dt,NT*dt])
legend({ ['N=' num2str(N(idx_min))], ...
               ['N=' num2str(N(idx_mid))], ...
               ['N=' num2str(N(idx_max))]}, 'Interpreter','latex', 'FontSize', fontSize, 'Location','northeast')
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', labelFontSize); 
ylabel('$z$ error [m]', 'Interpreter','latex', 'FontSize', labelFontSize)
    set(gca,'FontSize', fontSize);

%Zoom-in box (Manually located)
axInset = axes('Position', [0.4, 0.4, 0.3, 0.2]);  % [x y width height]
box on;
plot(tvec,z_mean_by_N(idx_min,:)+zeq-z_sim_trajectories(idx_min,:), 'color','m', 'LineWidth', lineWidth);hold on;grid on
plot(tvec,z_mean_by_N(idx_mid,:)+zeq-z_sim_trajectories(idx_mid,:), 'color',color2, 'LineWidth', lineWidth)
plot(tvec,z_mean_by_N(idx_max,:)+zeq-z_sim_trajectories(idx_max,:), 'color',color5, 'LineWidth', lineWidth)
xlim([1, 1.3]);   % zoomed x-range
ylim([-0.5e-4, 0.5e-4]);   % zoomed y-range
set(axInset, 'FontSize', fontSize);

% Draw connecting lines (normalized coordinates)
annotation('line', [0.3 0.162], [0.53 0.62], 'LineStyle', '--');
annotation('line', [0.3  0.162], [0.73 0.62], 'LineStyle', '--');

set(gcf, 'Units', 'centimeters', 'Position', [0 0 41 20]) % or [left bottom width height]
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [41 20])
set(gcf, 'PaperPositionMode', 'manual')
print(gcf, 'Figures/RMSE/mid_R_osc_zerror', '-dpdf', '-vector', '-fillpage');

savefig(folder + '/zerror_mid_R_osc_selectedN.fig')
% figure(3); clf
% plot(z_mean_by_N(idx_min,:)+zeq - z_sim_trajectories(idx_min,:), 'b--'); hold on
% plot(z_mean_by_N(idx_mid,:)+zeq - z_sim_trajectories(idx_mid,:), 'g--');
% plot(z_mean_by_N(idx_max,:)+zeq - z_sim_trajectories(idx_max,:), 'k--');
% yline(0, 'm--')
% legend({['N=' num2str(N(idx_min))], ...
%         ['N=' num2str(N(idx_mid))], ...
%         ['N=' num2str(N(idx_max))], 'zero error'})
% xlabel('Time Step'); ylabel('Estimation Error (z_{est} - z_{sim})')
% title('Estimation Error Trajectories')
% grid on
% savefig(folder + '/zs_error_selectedN.fig')

% figure(2);clf
% plot(z_sim_trajectories(7,:),'r-',LineWidth=2);hold on
% yline(zeq)
% legends=[];
% for i=1:size(z_mean_by_N,1)
%     plot(1:NT,z_mean_by_N(i,:)+zeq)
%     legends=[legends,'N = '+string(N(i))];
% end
% legend(["sim","zeq",legends])
% 
% savefig(folder + '/zs_all.fig')

%%
% figure(3);clf
% plot(mean(z_sim_trajectories(3,:),1)); hold on
% plot(z_mean_by_N(1,:)+zeq)
% plot(z_mean_by_N(7,:)+zeq)
% legend(["sim","N=6", "N=18"])
% savefig(folder + '/zs_reduced.fig')













%% 
save

%%

function [sim,est]=runSystem(N,dt,NT, X0, MHE_x0,zeq)

        params = parameters;

    
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
    
    t = 2;
    
    N_MHE = N;
    N_MPC = 10;
    %dt = 0.003;
    %NT = 350;
    tvec = 0:1:NT-1;
    
    %MHE tuning
    alpha = 0.9;
    noise_std = 0.1 * 1e-3; %0.1 mT
    R_MHE = inv(noise_std^2 * eye(nMeasurements));  
    Q_MHE=1e6 *diag([1,1,1,1,1,5,5,5,5,5]); 
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
            useAdvancedControl = true;
        end
        
        %disp(string(k) + ", Running with advanced control: " + ...
               % string(useAdvancedControl))
    
    
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
            
        %U = U + 10*sin(0.5*k).*[1.5;-1.5;-1.5;1.5];
        %U = zeros(4,1);

        %[~, X] = ode15s(@(t, x) f(x, U, params), tspan, X_sim(:,k));
        %X_sim(:, k+1) = X(end, :)'; 
        X_sim(:, k+1) = RK4Step(@f, X_sim(:,k), U, dt, params);
        U_sim(:, k) = U; 
        newU = U_sim(:, k); %For MHE input
    
        noise = noise_std * randn([nMeasurements, 1]);
        %noise=0;
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
        NEES_traj(k) = error' / mhe.currentP * error;
    
        %profile viewer
        elapsed = toc(t_start);
    end
    
    sim = X_sim;
    est = MHE_est;
end

function gamma = gamma_f(k,fade_period)
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end