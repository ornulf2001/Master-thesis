clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
import casadi.*

%Run msd simulation
run("msd_sim.m")


% Defining misc.
N_MHE = 8;
nStates = size(Ac,1);
nControls = size(Bc,2);
nMeasurements=size(C,1);
z0_block=zeros(nStates,1);

%Decision variables
X = SX.sym('X', nStates, N_MHE+1); %states
W = SX.sym('W', nStates, N_MHE); %Process noise
V = SX.sym('V', nMeasurements, N_MHE+1); %Measurement noise
z = [reshape(X, nStates*(N_MHE+1), 1);reshape(W, nStates*(N_MHE), 1);reshape(V, nMeasurements*(N_MHE+1), 1)];
P = zeros( nStates+nMeasurements+nControls, 1+(N_MHE+1)+N_MHE); 

 
R=0.001;
Q=diag([1,1]);
M=diag([1,1]);


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

xsol=zeros(nStates,size(Y_noisy,2)-(N_MHE+1));
iter_reg=0;
zs=[];

%Buffer and run regular MHE
%run("init_mhe.m")

%run("run_mhe.m")




%%%%%%%%%%%%%% Class version of MHE %%%%%%%%%%%%%%%%%%

mhe=MHEclass(N_MHE,Ac,Bc,C,Q,R,M,z0_block,x0_sim,dt);
iter_class=0;

%Buffer first horizon of measurements and control inputs
while ~mhe.isReadyToRun
    newY=Y_noisy(mhe.yBufferCount);
    newU=U_list(mhe.uBufferCount);
    mhe=mhe.bufferInitialData(newY, newU);
end


% Run class version of MHE
xsol2=zeros(nStates,size(Y_noisy,2)-(N_MHE+1));
 for k=1:size(Y_noisy,2)-(N_MHE+1)
     newY=Y_noisy(:,N_MHE+1+k);
     newU=U_list(N_MHE+k,:);
     mhe=mhe.runMHE(newY,newU);
     xsol2(:,k)=mhe.xCurrent;
     iter_class=iter_class+1;
     
 end
 

run("plotting.m")


