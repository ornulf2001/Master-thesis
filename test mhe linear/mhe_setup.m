clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")
import casadi.*
run("init_variables.m")

N_MHE = 4;

%Dynamical model
[A,B,z0_block]=f_func_lin(z_eq,params,nStates);




%Measurement model

%[C,D] = h_func_lin(z0,params);
C=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0];
nMeasurements=size(C,1);

%Decision variables
X = SX.sym('X', nStates, N_MHE+1); %Decision variables (states)
U = SX.sym('U', nControls,N_MHE); %Decision variables (controls)
W = SX.sym('W', nStates, N_MHE);
V = SX.sym('V', nMeasurements, N_MHE+1);

P = SX.sym('P', nStates+nMeasurements, 1+(N_MHE+1)); 
% %^ Parameters. Include the prior estimate x_(k-1) and N_MHE+1 y
% measurements (x,z,theta) NB Check this dimension!!! idk if we only store
% one measurement each step or all N+1 measurements while changing only one
% each time. I think we have first column for xprior from last step, then
% N+1 columns for the measurements of the whole horizon. Then each
% iteration we remove the oldest measurements and shift them, then append
% the newest measurement column at the end. This is the main complication
% of looking "backwards" in time in MHE, the rest of the code logic is the same as
% if looking "forward" in time I think. That is, make sure to update P
% correctly and remove/append the correct columns each time step. 


x_prior=P(1:nStates,1);
x_prior=zeros(size(x_prior));


for k = 1:(N_MHE+1)
    V(:,k) = P(nStates+1:end, k+1) - (C*X(:,k)); %+D*U(:,k)); % v_k = y_meas,k - h(x_k)
end

z=[reshape(X, nStates*(N_MHE+1), 1);reshape(U, nControls*(N_MHE), 1);reshape(W, nStates*(N_MHE), 1);reshape(V, nMeasurements*(N_MHE+1), 1)];
%^ z=[x0,z0,...thetadot0,...,xN,zN,...,thetadotN,   ux0,uz0,...,uxN-1,uzN-1 
% ,  w0,...wN-1,   v0,...,vN]
% 
meas_cov = diag([0.01^2, 0.01^2, deg2rad(3)^2]);            %noise covariance on measurements
proc_cov=diag([0.05^2, 0.05^2, deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %noise covariance on control input
arrival_cov=diag([0.05^2, 0.05^2, deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %Arrival cost covariance
L_y = chol(meas_cov,'lower'); R = L_y \ eye(size(L_y));     %These give R and Q as weight matrices for y and u resp. V=inv(sqrt(meas_cov)) W=...
L_x = chol(proc_cov, 'lower'); Q = L_x \eye(size(L_x));  
L_M = chol(arrival_cov, 'lower'); M = L_M\eye(size(L_M));

% Cost function quadratic terms
cost_X_block = blkdiag(M,zeros(nStates*N_MHE));
cost_U_block = kron(zeros(nControls),eye(N_MHE));
cost_Q_block = kron(Q,eye(N_MHE));
cost_R_block= kron(R,eye(N_MHE+1));
G=blkdiag(cost_X_block,cost_U_block,cost_Q_block,cost_R_block);

% Cost function linear terms
g_X=[-2*M*x_prior;zeros(nStates*N_MHE,1)]; %Only linear term is -2*M*x_prior*x(0)
g_U = zeros(nControls * N_MHE, 1); % No linear terms for U, W, or V
g_W = zeros(nStates * N_MHE, 1);
g_V = zeros(nMeasurements * (N_MHE + 1), 1);
g = [g_X; g_U; g_W; g_V];
g_num=full(g);
%Construct objective function
obj = z'*G*z + g'*z;

% The following is the implementation of the weird A_eq matrix and b_eq
% vector to be used in the QP equality constraint A_eq * z = b_eq
rows_Aeq= N_MHE*nStates;
cols_Aeq= (N_MHE+1)*nStates + N_MHE*nControls + N_MHE*nStates + (N_MHE+1)*nMeasurements;
Aeq=zeros(rows_Aeq,cols_Aeq);
beq = zeros(N_MHE * nStates, 1);
for k=0:N_MHE-1
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*k+1                                        :   nStates*(k+1))=-A; %Placeholder identities for A_k on blockdiagonal
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(k+1)+1                                    :   nStates*(k+1+1))=eye(nStates); %Identities on super-blockdiagonal
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(N_MHE+1)+nControls*k+1                    :   nStates*(N_MHE+1)+nControls*(k+1)) = -B;
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(N_MHE+1)+nControls*(N_MHE)+nStates*k+1    :   nStates*(N_MHE+1)+nControls*(N_MHE)+nStates*(k+1)) = -eye(nStates);
    
    beq(nStates * k + 1 : nStates * (k + 1)) = z0_block;
end




full_Aeq=full(Aeq); 
disp("yay")


%%%%%%%%%%%%%%%%% Solver settings and setup %%%%%%%%%%%%%%%%%%%%%%

% Define bounds for decision variables
lbx = -inf(size(z)); % Default to no lower bounds
ubx = inf(size(z));  % Default to no upper bounds

% Apply bounds to states (e.g., state limits)
for k = 0:N_MHE
    lbx(nStates*k+1:nStates*(k+1)) = [-2; -2; -pi/2; -10; -10; -pi]; % Example lower bounds
    ubx(nStates*k+1:nStates*(k+1)) = [ 2;  2;  pi/2;  10;  10;  pi]; % Example upper bounds
end

% Apply bounds to controls (e.g., actuator limits)
for k = 0:N_MHE-1
    lbx((N_MHE+1)*nStates + nControls*k+1 : (N_MHE+1)*nStates + nControls*(k+1)) = [-0.6; -0.6]; % Lower bounds
    ubx((N_MHE+1)*nStates + nControls*k+1 : (N_MHE+1)*nStates + nControls*(k+1)) = [ 0.6;  0.6]; % Upper bounds
end


opts = struct;
opts.qpsol = 'qpoases'; % Use qpOASES as the solver

% qpOASES-specific options
opts.qpsol_options.print_level = 'none';       % Suppress solver output (verbosity)
opts.qpsol_options.tol_stat = 1e-6;            % Stationarity tolerance
opts.qpsol_options.tol_eq   = 1e-6;            % Equality constraint tolerance
opts.qpsol_options.tol_ineq = 1e-6;            % Inequality constraint tolerance
opts.qpsol_options.tol_comp = 1e-6;            % Complementarity tolerance
opts.qpsol_options.max_iter = 100;             % Maximum number of iterations
%solver = qpsol('solver', 'qpoases', struct('x', z, 'f', obj, 'A', Aeq, 'lba', b_eq, 'uba', b_eq));

nWSR = 1000;
tic
[xOpt, fval, exitflag, iterations,lambda,auxOutput] = qpOASES(2*G, g, Aeq, -100000000*ones(length(z),1), 100000000*ones(length(z),1), beq, beq,nWSR);
toc