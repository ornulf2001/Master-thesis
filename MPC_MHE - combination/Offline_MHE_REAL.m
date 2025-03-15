clc,clear

addpath(genpath("NewModel_2D"))
addpath(genpath("simulator"))

%Dynamics
run("parameters.m")

index = @(A, i) A(i);
fz = @(z) index(f([0, z, zeros(1, 4)]', [0, 0]', params), 5);
zeq = fzero(fz, 0.1);
Xeq = [0, zeq, zeros(1, 4)]';
Ueq = [0, 0]';
[Ac, Bc, C] = linearizeModel(@f, @h, Xeq, Ueq, params);

%C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
xRef = [0; 0; 0; 0; 0; 0];
X0=[0;0.0005;0;0;0;0];
NT=500;

N_MHE=15;
N_MPC=20;
dt=0.002;

alpha=0.9;
load("noise_cov_matrix_real.mat") %R: Covariance of noise measurements "data_no_control.m"
%noise_std=0.1*1e-3; %mT
%R_MHE=inv(noise_std^2*eye(nMeasurements));         
R_MHE=inv(R);
Q_MHE=5e10*diag([1,1,1,1,1,1]);
M_MHE=5e8*diag([1,1,1,3,3,3]);
P0 = inv(M_MHE);

%load("data_no_control.mat")
load("data_control_and_perturbation.mat")

Y_noisy=[data.y.bx';data.y.bz']*1e-3;
Ip=data.u.Ix_plus';
In=data.u.Ix_minus';
U_list = [Ip-In;Ip+In];


MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,1e-5*Q_MHE,1e-5*R_MHE,1e-5*M_MHE,X0,P0,dt,MHE_options);

iter_class=0;
for k=1:size(Y_noisy,2)/4
    k
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