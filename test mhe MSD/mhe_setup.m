clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
import casadi.*

%Run msd simulation
run("msd_sim.m")


% Defining misc.
N_MHE = 4;
nStates = size(Ac,1);
nControls = size(Bc,2);
nMeasurements=size(C,1);
z0_block=zeros(nStates,1);

%Decision variables
X = SX.sym('X', nStates, N_MHE+1); %states
W = SX.sym('W', nStates, N_MHE); %Process noise
V = SX.sym('V', nMeasurements, N_MHE+1); %Measurement noise
z = [reshape(X, nStates*(N_MHE+1), 1);reshape(W, nStates*(N_MHE), 1);reshape(V, nMeasurements*(N_MHE+1), 1)];
%^ z=[x0,z0,...thetadot0,...,xN,zN,...,thetadotN,   ux0,uz0,...,uxN-1,uzN-1 
% ,  w0,...wN-1,   v0,...,vN]


P = zeros( nStates+nMeasurements+nControls, 1+(N_MHE+1)+N_MHE); 
% %^ Parameters. Include the prior estimate x_(k-1) and N_MHE+1 y
% measurements (x,z,theta)_k and the previous control input u_(k-1).


 

%meas_cov = diag(1); %, 0.01^2, deg2rad(3)^2]);            %noise covariance on measurements
%proc_cov=diag([0.3,0.3]);%, deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %noise covariance on control input
%arrival_cov=diag([0.1, 0.1]);% deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %Arrival cost covariance
%L_y = chol(meas_cov,'lower'); R = L_y \ eye(size(L_y));     %These give R and Q as weight matrices for y and u resp. V=inv(sqrt(meas_cov)) W=...
%L_x = chol(proc_cov, 'lower'); Q = L_x \eye(size(L_x));  
%L_M = chol(arrival_cov, 'lower'); M = L_M\eye(size(L_M));

R=0.003;
Q=diag([0.3,3]);
M=diag([0.05,0.03]);
%^ To surpress noise: Lower R

% Discretization %
dt=0.01;
A=expm(Ac*dt);
B=inv(Ac) * (A - eye(size(Ac))) * Bc;



Q=dt*Q;
R=dt*R;
M=dt*M;

% Cost function quadratic terms
cost_X_block = blkdiag(M,zeros(nStates*N_MHE));
cost_Q_block = kron(Q,eye(N_MHE));
cost_R_block= kron(R,eye(N_MHE+1));
G=blkdiag(cost_X_block,cost_Q_block,cost_R_block);

% Cost function linear terms
g_X=[-2*M*P(1:nStates,1);zeros(nStates*N_MHE,1)]; %Only linear term is -2*M*x_prior*x(0)
g_W = zeros(nStates * N_MHE, 1);
g_V = zeros(nMeasurements * (N_MHE + 1), 1);
g = [g_X; g_W; g_V];
g_num=full(g);

%Construct objective function
obj = z'*G*z + g'*z;

% The following is the implementation of the weird A_eq matrix and b_eq
% vector to be used in the QP equality constraint A_eq * z = b_eq
rows_Aeq= N_MHE*nStates+(N_MHE+1)*nMeasurements;
cols_Aeq= (N_MHE+1)*nStates + N_MHE*nStates + (N_MHE+1)*nMeasurements;
Aeq=zeros(rows_Aeq,cols_Aeq);
beq = zeros(rows_Aeq, 1);
for k = 0:N_MHE-1
    Aeq(nStates*k+1 : nStates*(k+1), nStates*k+1 : nStates*(k+1)) = -A;
    Aeq(nStates*k+1 : nStates*(k+1), nStates*(k+1)+1 : nStates*(k+2)) = eye(nStates);
    Aeq(nStates*k+1 : nStates*(k+1), nStates*(N_MHE+1) + nStates*k + 1 :nStates*(N_MHE+1) + nStates*(k+1)) = -eye(nStates);
end
for k=0:N_MHE
    Aeq(nStates*N_MHE+nMeasurements*k+1 : nStates*N_MHE+nMeasurements*(k+1), nStates*k+1 : nStates*(k+1))=C;
    Aeq(nStates*N_MHE+nMeasurements*k+1 : nStates*N_MHE+nMeasurements*(k+1),  nStates*(N_MHE+1)+nStates*(N_MHE)+nMeasurements*k+1 : nStates*(N_MHE+1)+nStates*(N_MHE)+nMeasurements*(k+1)) = eye(nMeasurements);
end

runtime=zeros(size(Y_noisy,2)-(N_MHE+1),1);
xsol=zeros(nStates,size(Y_noisy,2)-(N_MHE+1));
%%%%%%%%%%%%%%%%% Solver settings and setup %%%%%%%%%%%%%%%%%%%%%%


run("init_mhe.m")


mhe=MHEclass(N_MHE,Ac,Bc,C,z0_block,x0_sim,dt);
%class init

while ~mhe.isReadyToRun
    newY=Y_noisy(mhe.yBufferCount);
    newU=U_list(mhe.uBufferCount);
    mhe=mhe.bufferInitialData(newY, newU);
end

%run("run_mhe.m")
%run("plotting.m")


