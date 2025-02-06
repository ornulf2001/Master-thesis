clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")
import casadi.*

%Run msd simulation
run("msd_sim.m")


% Defining misc
N_MHE = 1;
nStates = size(A,1);
nControls = size(B,2);
nMeasurements=size(C,1);
z0_block=zeros(nStates,1);

%Decision variables
X = SX.sym('X', nStates, N_MHE+1); %states
W = SX.sym('W', nStates, N_MHE); %Process noise
V = SX.sym('V', nMeasurements, N_MHE+1); %Measurement noise

P = zeros( nStates+nMeasurements+nControls, 1+(N_MHE+1)+N_MHE); 
P(1:nStates,1)=x0_sim;
% %^ Parameters. Include the prior estimate x_(k-1) and N_MHE+1 y
% measurements (x,z,theta)_k and the previous control input u_(k-1).

z=[reshape(X, nStates*(N_MHE+1), 1);reshape(W, nStates*(N_MHE), 1);reshape(V, nMeasurements*(N_MHE+1), 1)];
%^ z=[x0,z0,...thetadot0,...,xN,zN,...,thetadotN,   ux0,uz0,...,uxN-1,uzN-1 
% ,  w0,...wN-1,   v0,...,vN]
 

%meas_cov = diag(1); %, 0.01^2, deg2rad(3)^2]);            %noise covariance on measurements
%proc_cov=diag([0.3,0.3]);%, deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %noise covariance on control input
%arrival_cov=diag([0.1, 0.1]);% deg2rad(2)^2, 0.05^2, 0.05^2, deg2rad(2)^2]); %Arrival cost covariance
%L_y = chol(meas_cov,'lower'); R = L_y \ eye(size(L_y));     %These give R and Q as weight matrices for y and u resp. V=inv(sqrt(meas_cov)) W=...
%L_x = chol(proc_cov, 'lower'); Q = L_x \eye(size(L_x));  
%L_M = chol(arrival_cov, 'lower'); M = L_M\eye(size(L_M));

R=0.01;
Q=diag([0.3,0.3]);
M=diag([0.1,0.1]);
%^ To surpress noise: Lower R

% Discretization %
dt=0.01;
A=expm(A*dt);
B=inv(A) * (expm(A*dt) - eye(size(A))) * B;
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


%%%%%%%%%%%%%%%%% Solver settings and setup %%%%%%%%%%%%%%%%%%%%%%

%Initialization, fill up the first horizon of measurements and control inputs
nWSR=1000; %Max iterations


for k=1:N_MHE+1 %Fill up first horizon of measurements before loop in P and then in beq
    P(nStates+1:nStates+nMeasurements, k+1)=Y_noisy(k);
    beq(nStates*N_MHE+nMeasurements*(k-1)+1 : nStates*N_MHE+nMeasurements*(k)) = P(nStates+1:nStates+nMeasurements, k+1);

end
for k=1:N_MHE %Fill up first horizon-1 of control inputs applied in P and then in beq
    P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+k)=U_list(k);
    beq(nStates*k-1 : nStates*(k)) = z0_block + B*P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,(N_MHE+1)+1+k: (N_MHE+1)+1+nControls*k);
end

xsol=[];
for k=1:size(Y_noisy,2)-(N_MHE+1)

[zOpt, fval, exitflag, iterations,lambda,auxOutput] = qpOASES(2*G, g, Aeq, -100000000*ones(length(z),1), 100000000*ones(length(z),1), beq, beq,nWSR); %Solve
xCurrent = zOpt(nStates*N_MHE + 1 : nStates*(N_MHE + 1));  % Extract and store X_N
xsol = [xsol, xCurrent];
P(1:nStates, 1) = xCurrent; % Update xprior for next iteration

%Shift measurement window
P(nStates+1:nStates+nMeasurements, 2:N_MHE+2)=[P(nStates+1:nStates+nMeasurements, 3:N_MHE+2),Y_noisy(N_MHE+1+k)]; 

%Shift control input window
P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+1:1+(N_MHE+1)+N_MHE)=[P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+1+1:1+(N_MHE+1)+N_MHE),U_list(k+N_MHE)]; 

%update the C*Xk + Vk = Ymeas,k constraint in beq with new measurement
beq(nStates*N_MHE+1 : nStates*N_MHE+nMeasurements*(N_MHE+1)) = P(nStates+1:nStates+nMeasurements, 2:N_MHE+2); 

%Update dynamics constraint with new control input in beq
beq(1:nStates*N_MHE) = [beq(1+nStates:nStates*N_MHE) ; z0_block + B*P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,(N_MHE+1)+1+nControls*N_MHE)]; 

% Update arrival cost with new xprior
g_X(1:nStates)=-2*M*P(1:nStates,1);
g = [g_X; g_W; g_V];
end



%%%%%%%% Plotting %%%%%%%%%
figure(1)
%plot(Y_noisy(1,N_MHE+1:end),"k") 
hold on
plot(xsol(1,:),"r"); hold on
plot(Y_sim(1,N_MHE+1:end),"b")
title("Estimated position. N: "+num2str(N_MHE)+", R: "+num2str(round(R,3))+", Q: "+mat2str(round(Q,3))+", M: "+mat2str(round(M,3)))
grid on
legend([ "est.","sim"])
xlabel("Time")
ylabel("x1")

figure(2)
plot(xsol(2,:)) 
hold on
plot(X_sim(2,N_MHE+1:end))
legend(["est.","sim"])
title("Estimated velcocity")
xlabel("Time")
ylabel("x2")
