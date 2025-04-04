clc,clear

addpath(genpath('3D model reduced order'))
%addpath(genpath("simulation_2"))

%Dynamics
%params=getParams();
params = parameters;

index = @(A, i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1 
zeq =  fzero(fz,0.1);
%[zeq, zEqInv, dzEq, dzEqInv] = computeSystemEquilibria(params,'fast');
%zeq=0.027;

xeq = [0,0,zeq,zeros(1,7)]';
xlp=xeq;
ueq = [0,0,0,0]';
ulp=ueq;
[Ac, Bc, C] = linearizeModel(@f, @h, xeq, ueq, params);

%f = @(x,u) maglevSystemDynamics(x,u,params,"fast");
%h = @(x,u) maglevSystemMeasurements(x,u,params,"fast");
% Define the point to linearize around
%xLp = [0,0,zEq(1),zeros(1,9)]'; % Linearizing around the equilibria
%uLp = zeros(length(params.solenoids.r),1);

%Linearization
%delta = 1e-6; % Step-size used in numerical linearization
%[Ac,Bc,C,~] = finiteDifferenceLinearization(f,h,xlp,ulp,delta);
% I = [1,2,3,4,5,7,8,9,10,11];
% Ac=Ac(I,I);
% Bc=Bc(I,:);
% C=C(:,I);
% xlp=xlp(I);


%C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
xRef = [0; 0; 0; 0; 0; 0];
X0=[0;0;zeq;0;0;0;0;0;0;0;];

N_MHE=10;
N_MPC=10;
dt=0.00845;



%MHE tuning
alpha=0.9;
noise_std=0.1*1e-3; %mT
R_MHE=inv(noise_std^2*eye(nMeasurements));  %Measurement noise weight = inv(measurement noise cov)      
%R_MHE=load("noise_cov.mat").R;
Q_MHE=1e6*diag([1e0,1e0,1e0,1e0,1e0,1e2,1e2,1e2,1e2,1e2]); 
    %Start out with low Q to trust measurements during start up, 
    %then increase Q after N_MHE+1. 
    %See below in loop
                                  

M_MHE = 1e5*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]); %Arrival cost weight initial guess (updates KF-style in loop)
P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling=1; %Scaling factor for better pos

%load("data_no_control.mat")
load("data_with_control_correct_current.mat")
%%
Y_noisy=[data.y.bx0';data.y.by0';data.y.bz0';data.y.bx1';data.y.by1';data.y.bz1';data.y.bx2';data.y.by2';data.y.bz2'];
%Ip=data.u.Ix_plus';
%In=data.u.Ix_minus';

U_list = [data.u.Ix_plus';data.u.Ix_minus';data.u.Iy_plus';data.u.Iy_minus'];

%%
MHE_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);
NT=size(Y_noisy,2)/2;
NIS_traj = zeros(NT-1,1);
Innovations_traj = zeros(nMeasurements,NT-1);
dof_NIS = mhe.nMeasurements;       % degrees of freedom (number of measurements)
alpha_NIS = 0.05;  % 95% confidence = 1 - alpha
lowerBound_NIS = chi2inv(alpha_NIS/2, dof_NIS);
upperBound_NIS = chi2inv(1 - alpha_NIS/2, dof_NIS);


iter_class=0;
for k=1:NT-1
    k
    % if k==mhe.N_MHE+2
    %     mhe.Q = 5e3*mhe.Q; %This relies on having enabled dynamic update of arrival cost in MHE. 
    %                        %That is, G must be updated with the new Q, which is done automatically when 
    %                        %updating M as well in arrival cost. If this is not done, we must update G here also.
    % end

    
    newY=Y_noisy(:,k+1);
    newU=U_list(:,k);
    mhe=mhe.runMHE(newY,newU);
    xsol2(:,k+1)=mhe.xCurrent;

    NIS_current=mhe.currentNIS;
    NIS_traj(k) = NIS_current;
    Innovations_traj(:,k)=mhe.currentInnovation;
    %a=[a,mhe.P];
    %run mpc(xk)
    %uk = zopt

    iter_class=iter_class+1;

end
%%
est_meas=C*(xsol2);

var=2;

figure(1)
clf
subplot(3,1,1)
plot(Y_noisy(var,1:NT-1)); hold on
plot(est_meas(var,1:NT-1))
legend(["meas","est"])
title("bx")

subplot(3,1,2)
plot(Y_noisy(var+1,1:NT-1)); hold on
plot(est_meas(var+1,1:NT-1))
legend(["meas","est"])
title("by")

subplot(3,1,3)
plot(Y_noisy(var+2,1:NT-1)); hold on
plot(est_meas(var+2,1:NT-1))
legend(["meas","est"])
title("bz")
% subplot(3,1,3)
% %plot(Y_noisy(3,1:ceil(size(Y_noisy,2)/4))); hold on
% plot(xsol2(3,:))
% legend("est")
% title("Theta")


%%

figure(5)
clf
subplot(3,1,1)
plot(xsol2(1,1:NT-1));
title("x")

subplot(3,1,2)
plot(xsol2(2,1:NT-1));
title("y")

subplot(3,1,3)
plot(xsol2(3,1:NT-1));
title("z")

figure(7)
clf
subplot(3,1,1)
plot(xsol2(6,1:NT-1));
title("xdot")

subplot(3,1,2)
plot(xsol2(7,1:NT-1));
title("ydot")

subplot(3,1,3)
plot(xsol2(8,1:NT-1));
title("zdot")
%%

figure(6)
clf
C1=C(1:9,1:5);
C2=C(1:9,6:10);
Cflip=[C2,C1];

est_meas_flip=Cflip*(xsol2+xlp);

plot(Y_noisy(3,1:NT-1));hold on
plot(est_meas_flip(3,1:NT-1));
%(1:9,1:5)=C()

%%
%figure(6)
%clf

%plot(Y_noisy(8,:))
figure(8);
clf
plot(NIS_traj, 'LineWidth', 1.5); hold on;
yline(lowerBound_NIS, '--r', 'LineWidth', 1.5);
yline(upperBound_NIS, '--r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('NIS');
title(['NIS trajectory with 95% Chi-square bounds (DoF = ' num2str(mhe.nMeasurements) ')']);
grid on;
legend('NIS', 'Lower 95% bound', 'Upper 95% bound');
function data = getParams()
    persistent params

    if isempty(params)
        parameters2;
    end
    data = params;
end