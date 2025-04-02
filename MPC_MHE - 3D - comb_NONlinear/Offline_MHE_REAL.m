clc,clear

addpath(genpath('3D model reduced order'))
addpath(genpath("simulator"))

%Dynamics
params=parameters;

index = @(A, i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1

zeq =  fzero(fz,0.1);

xeq = [0,0,zeq,zeros(1,7)]';
xlp=xeq;
ueq = [0,0,0,0]';
ulp=ueq;
[Ac, Bc, C] = linearizeModel(@f, @h, xeq, ueq, params);

%C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
xRef = [0; 0; 0; 0; 0; 0];
X0=[0.001;0.001;zeq;0;0;0;0;0;0;0;];
NT=500;

N_MHE=15;
N_MPC=20;
dt=0.002;

alpha=0.9;
%load("noise_cov_matrix_real.mat") %R: Covariance of noise measurements "data_no_control.m"
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));         
%R_MHE=inv(R);
Q_MHE=10e6*diag([100,100,2,100,100,100,100,1,100,100]); 
M_MHE = 1e2*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]);
P0 = inv(M_MHE);

%load("data_no_control.mat")
load("data_control_and_perturbation.mat")

Y_noisy=[data.y.bx';data.y.bz']*1e-3;
Ip=data.u.Ix_plus';
In=data.u.Ix_minus';
U_list = [Ip-In;Ip+In];


MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,xlp,P0,dt,MHE_options);

iter_class=0;
for k=1:size(Y_noisy,2)/2
    k
    if k==mhe.N_MHE+2
        mhe.Q = 5e3*mhe.Q; %This relies on having enabled dynamic update of arrival cost in MHE. 
                           %That is, G must be updated with the new Q, which is done automatically when 
                           %updating M as well in arrival cost. If this is not done, we must update G here also.
    end

    newY=Y_noisy(:,k);
    newU=U_list(:,k);
    mhe=mhe.runMHE(newY,newU);
    xsol2(:,k)=mhe.xCurrent;
    %a=[a,mhe.P];
    %run mpc(xk)
    %uk = zopt

    iter_class=iter_class+1;

end

est_meas=C*xsol2;



figure(1)
subplot(2,1,1)
plot(Y_noisy(1,1:ceil(size(Y_noisy,2)/4))); hold on
plot(est_meas(1,:))
legend(["meas","est"])
title("bx")

subplot(2,1,2)
plot(Y_noisy(2,1:ceil(size(Y_noisy,2)/4))-mean(Y_noisy(2,1:ceil(size(Y_noisy,2)/4)))); hold on
plot(est_meas(2,:)-mean(est_meas(2,:)))
legend(["meas","est"])
title("bz")
% subplot(3,1,3)
% %plot(Y_noisy(3,1:ceil(size(Y_noisy,2)/4))); hold on
% plot(xsol2(3,:))
% legend("est")
% title("Theta")