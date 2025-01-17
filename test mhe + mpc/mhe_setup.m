import casadi.*


%We synthesize now offline based on xx with range-bearing measurement model
meas_cov = diag([0.2^2, deg2rad(3)^2]);
numIter=length(xx(1,:));
r=zeros(numIter-1,1);
alpha=zeros(numIter-1,1);


for k=1:numIter-1
    r(k) = sqrt(xx(1,k)^2+xx(2,k)^2) + sqrt(meas_cov(1,1))*randn(1,1);
    alpha(k)= atan2(xx(2,k),xx(1,k)) + sqrt(meas_cov(2,2))*randn(1,1);
end
y_measurements = [r, alpha];



T = 0.2;
N_MHE = 4; %Now we do recursive receeding horizon %size(y_measurements,1)-1; % Estimation horizon, here horizon is full measurement array length run once

vMax = 0.6; vMin=-vMax;
omegaMax=pi/4; omegaMin= -omegaMax;

%Dynamical model
x=SX.sym('x'); y=SX.sym('s'); theta = SX.sym('theta');
states = [x;y;theta]; nStates=length(states);
v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; nControls=length(controls);
rhs = [v * cos(theta); v * sin(theta); omega];
f = Function('f', {states,controls}, {rhs}); %nonlinear mapping function f(x,u)

%Measurement model
r=SX.sym('r'); alpha = SX.sym('alpha'); % Range and bearing
measurement_rhs = [sqrt(x^2+y^2);atan2(y,x)];
h = Function('h', {states},{measurement_rhs});

%Decision variables
U = SX.sym('U', nControls,N_MHE); %Decision variables (controls)
X = SX.sym('X', nStates, N_MHE+1); %Decision variables (states)
P = SX.sym('P', 2, (N_MHE+1)+N_MHE); 
%^ Parameters. both y and u have 2 elements (x,y and v,omega) so we can collect 
%all measurements in a P matrix with 2 rows. First N_MHE+1 colums are for the
%y measurements, and the remaining N_MHE columns are for the u
%measurements. We want to estimate states aswell as what input was applied.
%These will be passed to the MHE.

L_y = chol(meas_cov,'lower'); V = L_y \ eye(size(L_y));  %These give V and W as weight matrices for y and u resp. V=inv(sqrt(meas_cov)) W=...
L_u = chol(con_cov, 'lower'); W = L_u \eye(size(L_u));  

%Construct objective function
obj = 0;
for k=1:N_MHE+1
    st=X(:,k); %Decision variable
    h_x=h(st); % h(x)
    y_tilde=P(:,k); % y measurements, first N_MHE+1 columns of P.
    obj=obj + (y_tilde-h_x)'*V*(y_tilde - h_x);
end
for k=1:N_MHE
    con_MHE= U(:,k);
    u_tilde = P(:,N_MHE+k);
    obj = obj + (u_tilde-con_MHE)' * W * (u_tilde-con_MHE);
end

%Add multiple shooting constraints to g
g= SX.sym('g', 3*N_MHE, 1);
for k=1:N_MHE
    st = X(:,k); con_MHE2 = U(:,k);
    st_next = X(:,k+1);
    f_value = f(st,con_MHE2);
    st_next_euler = st + (T*f_value);
    g(3*k-2:3*k) = (st_next - st_next_euler);
end
    
% MHE solver setup 
OPT_variables = [reshape(X, 3*(N_MHE+1), 1); reshape(U, 2*N_MHE, 1)];
nlp_mhe = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct();
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 0; %0, 3
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_mhe, opts);

args = struct();
args.lbg(1:3*(N_MHE)) = 0; % equality constraints
args.ubg(1:3*(N_MHE)) = 0; % equality constraints

args.lbx(1:3:3*(N_MHE+1), 1) = -2; % state x lower bound
args.ubx(1:3:3*(N_MHE+1), 1) = 2;  % state x upper bound
args.lbx(2:3:3*(N_MHE+1), 1) = -2; % state y lower bound
args.ubx(2:3:3*(N_MHE+1), 1) = 2;  % state y upper bound
args.lbx(3:3:3*(N_MHE+1), 1) = -pi/2; % state theta lower bound
args.ubx(3:3:3*(N_MHE+1), 1) = pi/2;  % state theta upper bound

args.lbx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE, 1) = vMin; % v lower bound
args.ubx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE, 1) = vMax; % v upper bound
args.lbx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE, 1) = omegaMin; % omega lower bound
args.ubx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE, 1) = omegaMax; % omega upper bound




