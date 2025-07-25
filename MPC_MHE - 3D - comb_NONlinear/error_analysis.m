clc; clear
addpath(genpath('3D model reduced order'))

N = 2:2:4;           % Range of N
dt = 0.003;          % Range of dt (can be vector later)

[NGrid, dtGrid] = meshgrid(N, dt);
[nDt, nN] = size(NGrid);
nSim = 2;            % Number of Monte Carlo simulations

NT = 500;
RMSE = zeros(nDt, nN, nSim);
sims = zeros(10,NT,numel(NGrid),nSim);
ests = zeros(10,NT,numel(NGrid),nSim);
for k = 1:nSim
    for i = 1:numel(NGrid)
        [row, col] = ind2sub(size(NGrid), i);
        currentN = NGrid(row, col);
        currentDt = dtGrid(row, col);

        disp("Simulation series " + k)
        disp("Running simulation with N = " + num2str(currentN) + " and dt = " + num2str(currentDt))

        tic
        [sim, est] = runSystem(currentN, currentDt,NT);
        sims(:,:,i,k) = sim;
        ests(:,:,i,k) = est;

        error = sim - est;
        RMSE(row, col, k) = sqrt(mean(error(:).^2));
        toc
    end
end

disp("Finished all simulations")

% Average across Monte Carlo runs
RMSE_mean = mean(RMSE, 3);
%%
folder="data_"+datestr(datetime("now"),"yyyymmdd_HHMMSS");
if ~exist(folder, 'dir')
    mkdir(folder);
end
save(folder+"/RMSEtotal","RMSE")
save(folder+"/sims","sims")
save(folder+"/ests","ests")
save(folder+"/RMSEmean","RMSE_mean")

%% plot
figure(1)
clf

try %Generally fails of length(dt)=1 -> RMSE is vector not grid
    surf(NGrid,dtGrid,RMSE_mean); grid on
    colorbar
    xlabel('Horizon Length N');
    ylabel('Time Step dt');
    zlabel('RMSE');
    title('Surface plot of RMSE vs. N and dt');

catch
    plot(N,RMSE_mean,"-r"); grid on; 
    xlabel('Horizon Length N');
    ylabel('RMSE');
    title('Plot of RMSE vs. N');
end
savefig('data/RMSE_vs_N.fig')











%%

function [sim,est]=runSystem(N,dt,NT)

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
    
    % Tuning
    xRef = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
    X0=[0;0;zeq;0;0;0;0;0;0;0;];
    
    N_MHE=N;
    N_MPC=20;
    
    
    % States: | x y z phi theta xdot ydot zdot phidot thetadot |
    
    alpha=0.7;
    noise_std=0.1*1e-3; %mT
    R_MHE=inv(noise_std^2*eye(nMeasurements));         
    Q_MHE=25e4*diag([100,100,100,100,100,100,100,100,100,100]); 
    
    M_MHE = 1e5*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
    P0 = inv(M_MHE);
    
    Q_MPC = diag([500 500 2000 10 10 1 1 10 1 1]);
    R_MPC = diag([0.2 0.2 0.2 0.2]);
    
    Q_LQR = diag([ ...
       1e1,1e1,1e1,1e1,1e1, ...
       1e1,1e1,1e5,1e1,1e1
       ]);
    R_LQR = 1e2*eye(4);
    
    
    % Bounds
    run("mpc_bounds.m")
    
    % Run
    MHE_options = optimoptions('quadprog','Display','none', 'Algorithm','interior-point-convex');
    mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,1e-5*Q_MHE,1e-5*R_MHE,1e-5*M_MHE,X0,xlp,P0,dt,MHE_options);
    
    MPC_options = optimoptions('quadprog','Display','none', 'Algorithm', 'interior-point-convex');
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
        try
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
                    [~, Uopt]=mpc.runMPC(X_sim(:,k));
                    Udiff=U_LQR - Uopt;
                    U = (1-gamma_f(k,fade_period))*U_LQR + gamma_f(k,fade_period)*Uopt;
                end
                
                [T, X] = ode45(@(t, x) f(x, U, params), tspan, X_sim(:,k));
                X_sim(:, k+1) = X(end, :)';
        
                U_sim(:,k) = U;
                %X_sim(:,k+1) = mpc.A*X_sim(:,k) + mpc.B*U_sim(:,k);
                newU=U_sim(:,k);
        
            else
                [~, Uopt]=mpc.runMPC(X_sim(:,k));
        
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
        
        
        catch
            continue
        end
    end
    
    sim=X_sim;
    est=MHE_est;
end

function gamma = gamma_f(k,fade_period)
    gamma = (k-fade_period(1))/(fade_period(end)-fade_period(1));
end