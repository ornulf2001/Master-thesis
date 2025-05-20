clc,clear
addpath(genpath('3D model reduced order_fixed'))

params = parameters;
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]', [0,0,0,0]', params), 8);  
zeq =  fzero(fz,0.1);
xeq = [0, 0, zeq, zeros(1,7)]'; xlp=xeq;
ueq = [0,0,0,0]';
N_MHE=10;
dt=0.003;
NT=350;
X0 = [0; 0; zeq; 0; 0; 0; 0; 0; 0; 0;];
MHE_x0 = X0-xlp;%zeros(nStates,1);

    
    
    
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
    
    N_MPC = 10;
    %dt = 0.003;
    %NT = 350;
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
    Q_MPC = diag([5 5 1000 1 1 1 1 200 1 1]);
    R_MPC = diag([0.002, 0.002, 0.002, 0.002]);
    
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
        k = iterCounter
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
        NEES_traj(k) = error' / mhe.currentP * error;
    
        %profile viewer
        elapsed = toc(t_start);
    end

    %%

    figure(1);clf
    plot(X_sim(3,:));hold on
    yline(zeq)