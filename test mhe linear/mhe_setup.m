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

%Construct objective function
obj = z'*G*z + g'*z;

% The following is the implementation of the weird A_eq matrix and b_eq
% vector to be used in the QP equality constraint A_eq * z = b_eq
rows_Aeq= N_MHE*nStates;
cols_Aeq= (N_MHE+1)*nStates + N_MHE*nControls + N_MHE*nStates + (N_MHE+1)*nMeasurements;
Aeq=zeros(rows_Aeq,cols_Aeq);
b_eq = zeros(N_MHE * nStates, 1);

% Compute z_0 offset contribution for z-dynamics


for k=0:N_MHE-1
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*k+1                                        :   nStates*(k+1))=-A; %Placeholder identities for A_k on blockdiagonal
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(k+1)+1                                    :   nStates*(k+1+1))=eye(nStates); %Identities on super-blockdiagonal
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(N_MHE+1)+nControls*k+1                    :   nStates*(N_MHE+1)+nControls*(k+1)) = -B;
    Aeq(nStates*k+1     :   nStates*(k+1), nStates*(N_MHE+1)+nControls*(N_MHE)+nStates*k+1    :   nStates*(N_MHE+1)+nControls*(N_MHE)+nStates*(k+1)) = -eye(nStates);
    
    b_eq(nStates * k + 1 : nStates * (k + 1)) = z0_block;
end




full_Aeq=full(Aeq); 
disp("yay")



% % MHE solver setup 
% OPT_variables = [reshape(X, 3*(N_MHE+1), 1); reshape(U, 2*N_MHE, 1)];
% nlp_mhe = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
% 
% opts = struct();
% opts.ipopt.max_iter = 2000;
% opts.ipopt.print_level = 0; %0, 3
% opts.print_time = 0;
% opts.ipopt.acceptable_tol = 1e-8;
% opts.ipopt.acceptable_obj_change_tol = 1e-6;
% 
% solver = nlpsol('solver', 'ipopt', nlp_mhe, opts);
% 
% args = struct();
% args.lbg(1:3*(N_MHE)) = 0; % equality constraints
% args.ubg(1:3*(N_MHE)) = 0; % equality constraints
% 
% args.lbx(1:3:3*(N_MHE+1), 1) = -2; % state x lower bound
% args.ubx(1:3:3*(N_MHE+1), 1) = 2;  % state x upper bound
% args.lbx(2:3:3*(N_MHE+1), 1) = -2; % state y lower bound
% args.ubx(2:3:3*(N_MHE+1), 1) = 2;  % state y upper bound
% args.lbx(3:3:3*(N_MHE+1), 1) = -pi/2; % state theta lower bound
% args.ubx(3:3:3*(N_MHE+1), 1) = pi/2;  % state theta upper bound
% 
% args.lbx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE, 1) = vMin; % v lower bound
% args.ubx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE, 1) = vMax; % v upper bound
% args.lbx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE, 1) = omegaMin; % omega lower bound
% args.ubx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE, 1) = omegaMax; % omega upper bound
% 
% 
% 
% 
function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end