clc,clear

addpath(genpath('3D model reduced order_fixed'))

%Dynamics
params = parameters2();

index = @(A, i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1 
zeq =  fzero(fz,0.03);
xeq = [0,0,zeq,zeros(1,7)]';
xlp=xeq;
ueq = [0,0,0,0]';
ulp=ueq;
[Ac, Bc, C] = linearizeModel(@f, @h, xeq, ueq, params);

%Tuning
%X0=[0.003;0.003;zeq+0.002;0;0;0;0;0;0;0;];
%X0=[0;0;0;0;0;0;0;0;0;0;];
X0=xlp;

alpha=1;
N_MHE=5;
dt=0.005;
A = expm(Ac * dt);
B = inv(Ac) * (A - eye(size(Ac))) * Bc;

%MHE tuning
noise_std=0.1*1e-3; %mT
%R_MHE=inv(noise_std^2*eye(size(C,1)));  %Measurement noise weight = inv(measurement noise cov)  
%R_MHE=1e2*eye(size(C,1));
%R_MHE=load("noise_cov.mat").R;
Q_MHE=1e5*diag([1e1,1e1,1e1,1e1,1e1,1e1,1e1,1e1,1e1,1e1]);                                   
M_MHE = 1e2*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]); %Arrival cost weight initial guess (updates KF-style in loop)
P0 = inv(M_MHE); % Arrival cost cov initial guess.
weightScaling =1;


%%
load("data_with_control_fixed_double_sample_rate.mat")
b0=data.sensorData{1};
b1=data.sensorData{2};
b2=data.sensorData{3};
%Y_noisy=[b0.bx';b0.by';b0.bz';b1.bx';b1.by';b1.bz';b2.bx';b2.by';b2.bz']*1e-3;

%U_list = [data.u.Ix_plus';data.u.Iy_plus';data.u.Ix_minus';data.u.Iy_minus'];
%U_list = [data.u.Ix_plus, data.u.Iy_plus, data.u.Ix_minus, data.u.Iy_minus]';
%U_list  = [data.u.Ix_plus,data.u.Ix_minus,data.u.Iy_plus,data.u.Iy_minus]';
% U_list = [data.u.Iy_minus';data.u.Ix_minus';data.u.Iy_plus';data.u.Ix_plus'];

%Y_noisy=load("Y_noisy_sim.mat").yNext_f;
%U_list =load("U_list_sim.mat").U_sim;
I = 100:700;

U_list = [data.u.Ix_plus(I), data.u.Iy_plus(I), data.u.Ix_minus(I), data.u.Iy_minus(I)]';
Y_noisy = 1e-3*[
    data.sensorData{1}.bx(I), data.sensorData{1}.by(I), data.sensorData{1}.bz(I),...
    data.sensorData{2}.bx(I), data.sensorData{2}.by(I), data.sensorData{2}.bz(I),...
    data.sensorData{3}.bx(I), data.sensorData{3}.by(I), data.sensorData{3}.bz(I)
]';
R_MHE=inv(cov(Y_noisy(:,400:end)'));
t=data.t(I);
t=t-t(1);

%% Run
MHE_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);
NT=ceil(size(Y_noisy,2));

xhat=X0;
yeq=h(xlp,ulp,params);

P_current=P0;
state_est=zeros(size(Ac,1),NT-1);
state_est(:,1)=xhat;
xsol2=zeros(size(Ac,1),NT-1);
vsol=zeros(mhe.nMeasurements,NT-1);
wsol=zeros(size(Ac,1),NT-1);
xsol2(:,1)=X0-xlp;
newY_f=Y_noisy(:,1);
for k=1:NT-1
    k
    
   %disp("P before shift") 
   %disp(mhe.P)

   

    newY=Y_noisy(:,k+1);
    newY_f=alpha*newY + (1-alpha)*newY_f; %EMA prefilter before MHE
    newY=newY_f; %For MHE input
    newU=U_list(:,k);

    mhe=mhe.runMHE(newY+C*xlp-yeq,newU);
    xsol2(:,k+1)=mhe.xCurrent;
    vsol(:,k+1)=mhe.vCurrent;
    wsol(:,k)=mhe.wCurrent;
    

    xhat_pred = A * xhat + B * newU;
    P_pred = A * P_current * A' + inv(Q_MHE);

    y_pred = yeq + C*(xhat_pred - xeq);
    innovation = newY - y_pred;
    S = C * P_pred * C' + inv(R_MHE);
    K = P_pred * C' / S;
    xhat = xhat_pred + K * innovation;
    P_current = (eye(size(P_pred)) - K * C) * P_pred * (eye(size(P_pred)) - K * C)' + K*inv(R_MHE)*K';
    state_est(:, k+1) = xhat;
   


end

est_meas = repmat(yeq, 1, NT) + C*(state_est - repmat(xeq, 1, NT));
est_meas2= repmat(yeq, 1, NT) + C*(xsol2 - repmat(xeq, 1, NT));

%%
%est_meas=C*(xsol2);

figure(1)
clf
plot(t,Y_noisy(1:6,:),'b-');hold on
plot(t,est_meas(1:6,:),'r-')
legend(["meas","est"])
title("KF")

figure(2)
clf
plot(t,Y_noisy(1:6,:),'b-');hold on
plot(t,est_meas2(1:6,:),'r-')
legend(["meas","est"])
title("MHE")
% subplot(3,1,1)
% plot(Y_noisy(4,1:NT-1)); hold on
% plot(est_meas(4,1:NT-1))
% legend(["meas","est"])
% title("b2x")
% 
% subplot(3,1,2)
% plot(Y_noisy(5,1:NT-1)); hold on
% plot(est_meas(5,1:NT-1))
% legend(["meas","est"])
% title("b2y")
% 
% subplot(3,1,3)
% plot(Y_noisy(6,1:NT-1)); hold on
% plot(est_meas(6,1:NT-1))
% legend(["meas","est"])
% title("b2z")
% % 
% 
% subplot(9,1,4)
% plot(Y_noisy(4,1:NT-1)); hold on
% plot(est_meas(4,1:NT-1))
% legend(["meas","est"])
% title("b1x")
% 
% subplot(9,1,5)
% plot(Y_noisy(5,1:NT-1)); hold on
% plot(est_meas(5,1:NT-1))
% legend(["meas","est"])
% title("b1y")
% 
% subplot(9,1,6)
% plot(Y_noisy(6,1:NT-1)); hold on
% plot(est_meas(6,1:NT-1))
% legend(["meas","est"])
% title("b1z")
% 
% 
% subplot(9,1,7)
% plot(Y_noisy(7,1:NT-1)); hold on
% plot(est_meas(7,1:NT-1))
% legend(["meas","est"])
% title("b2x")
% 
% subplot(9,1,8)
% plot(Y_noisy(8,1:NT-1)); hold on
% plot(est_meas(8,1:NT-1))
% legend(["meas","est"])
% title("b2y")
% 
% subplot(9,1,9)
% plot(Y_noisy(9,1:NT-1)); hold on
% plot(est_meas(9,1:NT-1))
% legend(["meas","est"])
% title("b2z")
% 
% 
% %%
% xsol2=xsol2+xlp;
% figure(5)
% clf
% subplot(3,1,1)
% plot(xsol2(1,1:NT-1));
% title("x")
% 
% subplot(3,1,2)
% plot(xsol2(2,1:NT-1));
% title("y")
% 
% subplot(3,1,3)
% plot(xsol2(3,1:NT-1));
% title("z")
% 
% figure(7)
% clf
% subplot(3,1,1)
% plot(xsol2(6,1:NT-1));
% title("xdot")
% 
% subplot(3,1,2)
% plot(xsol2(7,1:NT-1));
% title("ydot")
% 
% subplot(3,1,3)
% plot(xsol2(8,1:NT-1));
% title("zdot")
% %%
% 
% % figure(6)
% % clf
% % C1=C(1:9,1:5);
% % C2=C(1:9,6:10);
% % Cflip=[C2,C1];
% % 
% % est_meas_flip=Cflip*(xsol2);
% % 
% % plot(Y_noisy(8,1:NT-1));hold on
% % plot(est_meas_flip(3,1:NT-1));
% %(1:9,1:5)=C()
% 
% %%
% % figure(6)
% % clf
% % plot(Y_noisy(8,:))
% figure(8);
% clf
% plot(NIS_traj, 'LineWidth', 1.5); hold on;
% yline(lowerBound_NIS, '--r', 'LineWidth', 1.5);
% yline(upperBound_NIS, '--r', 'LineWidth', 1.5);
% xlabel('Time');
% ylabel('NIS');
% title(['NIS trajectory with 95% Chi-square bounds (DoF = ' num2str(mhe.nMeasurements) ')']);
% grid on;
% legend('NIS', 'Lower 95% bound', 'Upper 95% bound');
% 
% %%
% 
% est_meas2=C*(state_est+xlp);
% 
% figure(9)
% clf
% subplot(3,1,1)
% plot(Y_noisy(4,1:NT-1)); hold on
% plot(est_meas2(4,1:NT-1))
% legend(["meas","est"])
% title("b0x")
% 
% subplot(3,1,2)
% plot(Y_noisy(5,1:NT-1)); hold on
% plot(est_meas2(5,1:NT-1))
% legend(["meas","est"])
% title("b0y")
% 
% subplot(3,1,3)
% plot(Y_noisy(6,1:NT-1)); hold on
% plot(est_meas2(6,1:NT-1))
% legend(["meas","est"])
% title("b0z")
% 
% 
% subplot(9,1,4)
% plot(Y_noisy(4,1:NT-1)); hold on
% plot(est_meas2(4,1:NT-1))
% legend(["meas","est"])
% title("b1x")
% 
% subplot(9,1,5)
% plot(Y_noisy(5,1:NT-1)); hold on
% plot(est_meas2(5,1:NT-1))
% legend(["meas","est"])
% title("b1y")
% 
% subplot(9,1,6)
% plot(Y_noisy(6,1:NT-1)); hold on
% plot(est_meas2(6,1:NT-1))
% legend(["meas","est"])
% title("b1z")
% 
% 
% subplot(9,1,7)
% plot(Y_noisy(7,1:NT-1)); hold on
% plot(est_meas2(7,1:NT-1))
% legend(["meas","est"])
% title("b2x")
% 
% subplot(9,1,8)
% plot(Y_noisy(8,1:NT-1)); hold on
% plot(est_meas2(8,1:NT-1))
% legend(["meas","est"])
% title("b2y")
% 
% subplot(9,1,9)
% plot(Y_noisy(9,1:NT-1)); hold on
% plot(est_meas2(9,1:NT-1))
% legend(["meas","est"])
% title("b2z")
%%

% figure(101);clf
% plot(state_est(3,:)+xlp(3));hold on
% plot(xsol2(3,:)+xlp(3))
% title("z")
% legend(["KF","MHE"])
% 
% figure(102);clf
% plot(state_est(8,:)+xlp(8));hold on
% plot(xsol2(8,:)+xlp(8))
% title("zdot")
% legend(["KF","MHE"])
% 
% figure(103);clf
% plot(state_est(2,:)+xlp(2));hold on
% plot(xsol2(2,:)+xlp(2))
% title("x")
% legend(["KF","MHE"])
% 
% figure(104);clf
% plot(state_est(7,:)+xlp(7));hold on
% plot(xsol2(7,:)+xlp(7))
% title("xdot")
% legend(["KF","MHE"])






%%
function data = getParams()
    persistent params

    if isempty(params)
        parameters2;
    end
    data = params;
end