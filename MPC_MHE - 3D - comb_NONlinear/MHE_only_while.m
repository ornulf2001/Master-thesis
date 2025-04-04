clc,clear
%addpath(genpath("3D model"))
addpath(genpath('3D model reduced order'))
addpath(genpath('../qpOASES/interfaces/matlab'))


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

% States: [ x y z phi theta xdot ydot zdot phidot thetadot ]
[Ac,Bc,C] = linearizeModel(@f,@h,xlp,ulp,params);


nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%X0=[0.001;0.001;zeq;0;0;0;0;0;0;0;];
X0 = zeros(nStates,1);
NT=200;

N_MHE=10;
dt=0.003;

alpha=0.9;
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));  %Measurement noise weight = inv(measurement noise cov)      
Q_MHE=10e6*diag([100,100,10,10,100,100,100,100,100,10]); 
    %Start out with low Q to trust measurements during start up, 
    %then increase Q after N_MHE+1. 
    %See below in loop
                                  

M_MHE = 1e2*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]); %Arrival cost weight initial guess (updates KF-style in loop)
P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling=1e-4; %Scaling factor for better posing of hessian


MHE_options = optimset('Display','on', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);

MHE_est = zeros(nStates, NT);
MHE_est(:,1)=mhe.x0; xEst = mhe.x0;
U_sim = zeros(nControls, NT-1);
yNext=zeros(nMeasurements,NT);  
%yNext(:,1)= C*(X0-xlp);
yNext_f=zeros(nMeasurements,NT);
%yNext_f(:,1)=%C*(X0-xlp);
NIS_traj = zeros(NT-1,1);
NEES_traj = zeros(NT-1,1);
Innovations_traj = zeros(nMeasurements,NT-1);
newY=yNext(:,1);
iterCounter = 1;

while iterCounter<NT
    disp(iterCounter)
    k=iterCounter;
    if iterCounter==mhe.N_MHE+2
        mhe.Q = 5e3*mhe.Q;
        mhe.G(mhe.nStates*(mhe.N_MHE+1)+1:mhe.nStates*(mhe.N_MHE+1)+mhe.nStates*mhe.N_MHE ,mhe.nStates*(mhe.N_MHE+1)+1:mhe.nStates*(mhe.N_MHE+1)+mhe.nStates*mhe.N_MHE ) = kron(mhe.weightScaling*mhe.Q,eye(mhe.N_MHE));
        %Here we increase Q after N_MHE+1 iterations when the MHE has calibrated. This seems to improve performance?
    end


    U_sim(:,k)=[0;0;0;0];
    newU=U_sim(:,k); %For MHE input

    noise=noise_std*randn([nMeasurements,1]);
    %yNext(:,k+1) = C*X_sim(:,k+1)-C*xlp+noise; %Subtract xlp to correct the frame of ref (?)
    yNext(:,k+1)=zeros(nMeasurements,1) + noise;
    yNext_f(:,k+1)=alpha*yNext(:,k+1) + (1-alpha)*yNext_f(:,k); %EMA prefilter before MHE
    newY=yNext(:,k+1); %For MHE input

    mhe=mhe.runMHE(newY,newU); %Estimate xk with MHE
    xEst=mhe.xCurrent; %xk^
    MHE_est(:,k+1)=xEst;
    iterCounter=iterCounter+1;
end

%%
mhe_meas=C*MHE_est;
figure(1)
clf
plot(yNext(1,:));hold on
plot(mhe_meas(1,:)); 
legend(["meas","est"])
