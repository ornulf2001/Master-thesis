clc,clear
% addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
% addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("NewModel_2D"))
% import casadi.*

%Dynamics
params=getParams();
index = @(A,i) A(i);
fz = @(z) index(f([0,z,zeros(1,4)]',[0,0]',params),5);

zeq =  fzero(fz,0.1);
Xeq = [0, zeq, zeros(1, 4)]';
Ueq = [0, 0]';
[Ac,Bc,C] = linearizeModel(@f,@h,Xeq,Ueq,params);
% z_eq=0.0365;
% Xeq = [0; z_eq; 0; 0; 0; 0];
% Ueq=[0;0];
% delta = 1e-6;
% f = maglevSystemDynamics2d(x,u,params);
% h = @(x,u) maglevSystemMeasurements2d(x,u,params);
% [Ac, Bc, ~, ~] = finiteDifferenceLinearization(f, h, Xeq, Ueq, delta);
% C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
X0=[0;0;0;0;0;0];
xRef = [0; 0.03; 0; 0; 0; 0];
NT=500;

N_MHE=15;
N_MPC=20;
dt=0.003;

R_MHE=5*dt*diag([1,1]);
Q_MHE=10*diag([1,1,1,1,1,1]);
M_MHE=0.1*diag([1,1,1,3,3,3]);
noise_cov=0.01;

Q_MPC = diag([5000 1000 5000 1 100 1]);
R_MPC = diag([0.1 0.01]);

%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','on', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, nStates, nControls, MPC_options, xRef, lbuRef, ubuRef);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT-1);
MHE_est = zeros(nStates, NT);
yNext=zeros(nMeasurements,NT);  
yNext(:,1)= C*X0;
newY=yNext(:,1);
xNext = X0;
MPC_Xopt(:, 1) = X0;
xCurrent = X0;
error=[];

tspan = [0, dt];

% Calculating the reference input for stabilizing in the reference point
uRef = mpc.computeReferenceInput();

init_plot(X0);

for k=1:NT
  
    if k<50
        U=[0;0];
        % MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;  
        [T, X] = ode15s(@(t, x) f(x, U, params), tspan, MPC_Xopt(:, k));
        
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        xNext = MPC_Xopt(:, k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);
        xCurrent = xNext;
        

    elseif k>=50 && k<60
        U=[1;3];
        % MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        [T, X] = ode15s(@(t, x) f(x, U, params), tspan, xEst);
       
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        xNext = MPC_Xopt(:, k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);
        xCurrent = xNext;
        MPC_Xtrue(:, k+1) = xCurrent;

    else
        
        [xTrue, Uopt]=mpc.runMPC(xEst);
        xCurrent = xTrue;
        
        [T, X] = ode15s(@(t, x) f(x, Uopt, params), tspan, xEst);
        
        MPC_Xopt(:, k+1) = X(end, :)';% + Xeq;
        xNext = MPC_Xopt(:, k+1);
        MPC_Uopt(:,k) = Uopt;% + uRef;
        newU = MPC_Uopt(:,k);
    end
    noise=noise_cov*dt*randn(size(C*xNext));
    yNext(:,k+1) = C*xNext + noise;
    newY=yNext(:,k+1);
    mhe=mhe.runMHE(newY,newU);
    xEst = mhe.xCurrent + Xeq; % This is the important aspect, changing coordinates for the estimates before feeding it to the simulation
    MHE_est(:,k)=xEst;
    error=[error,xNext-xEst];
    update_plot(xEst)
end




% figure(4)
% plot(MHE_est(1,1:end)); hold on
% plot(MPC_Xopt(1, 1:end));
% plot(yNext(1,:))
% title("X")
% legend('est', 'xopt',"meas")
% 
% figure(5)
% plot(MHE_est(2,1:end)); hold on
% plot(MPC_Xopt(2, 1:end));
% title("Z")
% legend('est', 'xopt')
% 
% figure(6)
% plot(MHE_est(3,1:end)); hold on
% plot(MPC_Xopt(3, 1:end));
% title("Theta")
% legend('est', 'xopt')
mpc.plotResults(MPC_Xopt,MPC_Uopt)



% Plotting Functions
function init_plot(x0)
    % Persistent variables for plotting
    MaglevPlotFig = getappdata(0, 'MaglevPlotFig');

    if isempty(MaglevPlotFig) || ~ishandle(MaglevPlotFig)
        % Create figure and axes
        MaglevPlotFig = figure('Name', 'Maglev Animation', 'NumberTitle', 'off');
        movegui(MaglevPlotFig, 'center');
        MaglevPlotAxes = axes('Parent', MaglevPlotFig);
        grid(MaglevPlotAxes, 'on'); hold(MaglevPlotAxes, 'on'); box(MaglevPlotAxes, 'on');
        daspect(MaglevPlotAxes, [1, 1, 1]);
        xlim(MaglevPlotAxes, [-0.1, 0.1]);
        ylim(MaglevPlotAxes, [-0.05, 0.3]);
        xlabel(MaglevPlotAxes, '$x$','Interpreter','latex','FontSize',14);
        ylabel(MaglevPlotAxes, '$z$','Interpreter','latex','FontSize',14);

        % Draw ground
        yline(MaglevPlotAxes, 0, 'k', 'LineWidth', 2);

        % Draw solenoids
        w = 2*0.0092;
        h = 0.0120;
        rectangle('Parent', MaglevPlotAxes, 'Position',[0.02-w/2,0,w,h], 'EdgeColor', 'k', 'FaceColor',[0.72, 0.45, 0.2], 'LineWidth',2);
        rectangle('Parent', MaglevPlotAxes, 'Position',[-0.02-w/2,0,w,h], 'EdgeColor', 'k', 'FaceColor',[0.72, 0.45, 0.2], 'LineWidth',2);

        % Draw maglev vehicle
        MaglevHandle = create_maglev(MaglevPlotAxes);
        set(MaglevHandle, 'Matrix', makehgtform('translate', [x0(1), x0(2), 0], 'zrotate', -x0(3)));

        % Store handles in application data
        setappdata(0, 'MaglevPlotFig', MaglevPlotFig);
        setappdata(0, 'MaglevHandle', MaglevHandle);
        setappdata(0, 'LastUpdateTime', 0);
    end
end

function update_plot(x)
    % Retrieve persistent data
    MaglevHandle = getappdata(0, 'MaglevHandle');

    % Update maglev position and orientation
    set(MaglevHandle, 'Matrix', makehgtform('translate', [x(1), x(2), 0], 'zrotate', x(3)));
    drawnow limitrate nocallbacks; % Efficient drawing
end

function MaglevHandle = create_maglev(ax)
    % Derived parameters
    w_body = 2*0.0250;
    h_body = 0.0050;
    r_center = 0.001;
    linewidth = 1;

    % Maglev transform object
    MaglevHandle = hgtransform('Parent', ax);

    % Body
    rectangle('Parent', MaglevHandle, 'Position', [-w_body/2, -h_body/2, w_body, h_body], ...
        'EdgeColor', 'k', 'FaceColor', [0.5, 0.5, 0.5], 'LineWidth', linewidth);

    % Center marker
    rectangle('Parent', MaglevHandle, 'Position', [-r_center/2, -r_center/2, r_center, r_center], ...
        'Curvature', [1, 1], 'EdgeColor', 'k', 'FaceColor', 'k', 'LineWidth', linewidth);
end