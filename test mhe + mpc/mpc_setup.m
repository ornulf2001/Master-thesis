
import casadi.*
%MPC implementation
T = 0.2; %[s]
N = 15; % Prediction horizon
robDiam = 0.1;

vMax = 0.6; vMin=-vMax;
omegaMax=pi/4; omegaMin= -omegaMax;

x=SX.sym('x'); y=SX.sym('s'); theta = SX.sym('theta');
states = [x;y;theta]; nStates=length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; nControls=length(controls);
rhs = [v * cos(theta); v * sin(theta); omega];%+ sqrt(W)*randn(3,1); %system rhs

f = Function('f', {states,controls}, {rhs}); %nonlinear mapping function f(x,u)
U = SX.sym('U', nControls,N); %Decision variables (controls)
P = SX.sym('P', nStates + nStates); %Parameters (which includes the initial state of the robot and the reference state)

X = SX.sym('X', nStates, N+1); %Vector that represents the states over the optimization problem.

obj = 0; %Objective function
g = SX.sym('g', 3+3*N+N+1, 1); %constraints
Q = zeros(3,3); Q(1,1) = 8; Q(2,2) = 8; Q(3,3) = 1; % State weigths
R = zeros(2,2); R(1,1) = 0.001; R(2,2) = 0.003; % Controls weights

con_cov=diag([0.05^2, deg2rad(2)^2]); %noise on control input

st = X(:,1); %initial state
g(1:3) = st-P(1:3); %add initial condition eq constraints

for k=1:N
    st = X(:,k); con = U(:,k);
    obj = obj + (st-P(4:6))' * Q * (st-P(4:6)) + con' * R * con; %Calculate objective function value
    stNext= X(:,k+1);
    fValue= f(st,con);
    stNextEuler = st + T*fValue;
    g(3*k+1:3*k+3) = stNext - stNextEuler; % add dynamical model eq constraints
end

%Add ineq constraint for collision avoidance with sphere obstacle

obsX=0.5; %meters
obsY=0.5; %meters
obsDiam= 0.3; %meters
for k=1:N+1 %box constraints due to the map margins
    g(3+3*N+k)= -((X(1,k)-obsX)^2 + (X(2,k)-obsY)^2) + (robDiam/2 + obsDiam/2)^2;
end

%Make optimization variables as a column vector
OptVariables = [reshape(X,nStates*(N+1),1); reshape(U,nControls*N,1)];
nlp_prob = struct('f', obj, 'x', OptVariables, 'g', g, 'p', P); % Organizing the optimization problem stuff in a struct

%Solver & options
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level=0; % 0, 3
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-3;
opts.ipopt.acceptable_obj_change_tol = 1e-3;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

%Optimization args
args = struct;
args.lbg(1:3*(N+1)) = 0; %equality constraints
args.ubg(1:3*(N+1)) = 0; %equality constraints

args.lbg(3*(N+1)+1 : 3*(N+1)+(N+1)) = -inf; %inequality constraints
args.ubg(3*(N+1)+1 : 3*(N+1)+(N+1)) = 0; %inequality constraints

args.lbx(1:3:3*(N+1),1) = -2; %state x lower bound
args.ubx(1:3:3*(N+1),1) =  2; %state x upper bound
args.lbx(2:3:3*(N+1),1) = -2; %state y lower bound
args.ubx(2:3:3*(N+1),1) =  2; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) =  inf; %state theta upper bound

args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = vMin; %v lower bound
args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = vMax; %v upper bound
args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omegaMin; %omega lower bound
args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omegaMax; %omega upper bound

% THE SIMULATION LOOP STARTS HERE
%-----------------------------------------------------------------------

t0 = 0;
x0 = [0;0;0]; % Initial conditon
xs = [1.2; 1.2; pi/4]; %Setpoint pose
xx(:,1) = x0; % xx contains the history of the states
u0 = zeros(N,2); % Two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % Initialization of the states decision variables
simTime= 20;

% Trigger MPC
mpcIter= 0;
maxIter=simTime/T;
xx1=zeros(N+1,length(states), maxIter);
u_cl=zeros(maxIter,2);
con=zeros(maxIter,2);
t=zeros(maxIter,1);
t(1)=t0;
