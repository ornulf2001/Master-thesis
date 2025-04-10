clc,clear

addpath(genpath('3D model reduced order_fixed'))
%addpath(genpath("simulation_2"))

%Dynamics
%params=getParams();
params = parameters();

index = @(A, i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1 
zeq =  fzero(fz,0.03);
% [zeq, zEqInv, dzEq, dzEqInv] = computeSystemEquilibria(params,'fast');
%zeq=0.027;

xeq = [0,0,zeq,zeros(1,7)]';
xlp=xeq;
ueq = [0,0,0,0]';
ulp=ueq;
%[Ac, Bc, C] = linearizeModel(@f, @h, xeq, ueq, params);
load("ABC_simple_model_reduced.mat")
% f = @(x,u) maglevSystemDynamics(x,u,params,"fast");
% h = @(x,u) maglevSystemMeasurements(x,u,params,"fast");
% Define the point to linearize around
%xLp = [0,0,zEq(1),zeros(1,9)]'; % Linearizing around the equilibria
%uLp = zeros(length(params.solenoids.r),1);


%Linearization
% delta = 1e-6; % Step-size used in numerical linearization
% [Ac,Bc,C,~] = finiteDifferenceLinearization(f,h,xlp,ulp,delta);
% I = [1,2,3,4,5,7,8,9,10,11];
% Ac=Ac(I,I);
% Bc=Bc(I,:);
% C=C(:,I);
% xlp=xlp(I);

Ac=A;Bc=B;

D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
xRef = [0; 0; 0; 0; 0; 0];
X0=[0;0;0;0;0;0;0;0;0;0;];

N_MHE=15;
dt=0.00845;

 A = expm(Ac * dt);
 B = inv(Ac) * (A - eye(size(Ac))) * Bc;
%load("data_no_control.mat")
load("data_with_control_correct_current.mat")

%MHE tuning
alpha=0.9;
%noise_std=0.1*1e-3; %mT
%R_MHE=inv(noise_std^2*eye(nMeasurements));  %Measurement noise weight = inv(measurement noise cov)      
%R_MHE=1e1*load("noise_cov.mat").R;
R_MHE=1e1*load("noise_cov.mat").R;
Q_MHE=1e3*diag([1e2,1e2,1e2,1e2,1e2,1e3,1e3,1e3,1e3,1e3]);                                   

M_MHE = 1e1*diag([5,5,5,0.005,005,0.002,0.002,0.002,0.0001,0.0001]); %Arrival cost weight initial guess (updates KF-style in loop)
P0 = inv(M_MHE); % Arrival cost cov initial guess.


%%
Y_noisy=[data.y.bx0';data.y.by0';data.y.bz0';data.y.bx1';data.y.by1';data.y.bz1';data.y.bx2';data.y.by2';data.y.bz2']*1e-3;
%R_MHE=inv(cov(Y_noisy'));

%U_list = [data.u.Ix_plus';data.u.Iy_plus';data.u.Ix_minus';data.u.Iy_minus'];
%U_list = [data.u.Ix_plus, data.u.Iy_plus, data.u.Ix_minus, data.u.Iy_minus]';
%U_list = [data.u.Iy_minus';data.u.Ix_minus';data.u.Iy_plus';data.u.Ix_plus'];

%% Run
%MHE_options = optimset('Display','off', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'active-set');
%mhe = MHEclass_KF_Update(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,weightScaling,X0,xlp,P0,dt,MHE_options);
NT=ceil(size(Y_noisy,2)/2);
% NIS_traj = zeros(NT-1,1);
% Innovations_traj = zeros(nMeasurements,NT-1);
% dof_NIS = mhe.nMeasurements;       % degrees of freedom (number of measurements)
% alpha_NIS = 0.05;  % 95% confidence = 1 - alpha
% lowerBound_NIS = chi2inv(alpha_NIS/2, dof_NIS);
% upperBound_NIS = chi2inv(1 - alpha_NIS/2, dof_NIS);
xhat=X0;
P_current=P0;

iter_class=0;
for k=1:NT-1
    
    
   

    
    %newY=Y_noisy(:,k+1);
    %newU=U_list(:,k);
    %mhe=mhe.runMHE(newY,newU);
    %xsol2(:,k+1)=mhe.xCurrent;

    u=U_list(:,k);
    xhat_pred = A * xhat + B * u;
    P_pred = A * P_current * A' + inv(Q_MHE);

    y = Y_noisy(:, k+1);
    innovation = y - C * xhat_pred;
    S = C * P_pred * C' + inv(R_MHE);
    K = P_pred * C' / S;
    xhat = xhat_pred + K * innovation;
    P_current = (eye(size(P_pred)) - K * C) * P_pred;
    state_est(:, k) = xhat;
    % NIS_current=mhe.currentNIS;
    % NIS_traj(k) = NIS_current;
    % Innovations_traj(:,k)=mhe.currentInnovation;

    iter_class=iter_class+1;

end



%%
% est_meas=C*(xsol2);
% Y_noisy=Y_noisy+C*xlp;
% 
% figure(1)
% clf
% subplot(9,1,1)
% plot(Y_noisy(1,1:NT-1)); hold on
% plot(est_meas(1,1:NT-1))
% legend(["meas","est"])
% title("b0x")
% 
% subplot(9,1,2)
% plot(Y_noisy(2,1:NT-1)); hold on
% plot(est_meas(2,1:NT-1))
% legend(["meas","est"])
% title("b0y")
% 
% subplot(9,1,3)
% plot(Y_noisy(3,1:NT-1)); hold on
% plot(est_meas(3,1:NT-1))
% legend(["meas","est"])
% title("b0z")
% 
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
est_meas2=C*(state_est);

figure(9)
clf
subplot(9,1,1)
plot(Y_noisy(1,1:NT-1)); hold on
plot(est_meas2(1,1:NT-1))
legend(["meas","est"])
title("b0x")

subplot(9,1,2)
plot(Y_noisy(2,1:NT-1)); hold on
plot(est_meas2(2,1:NT-1))
legend(["meas","est"])
title("b0y")

subplot(9,1,3)
plot(Y_noisy(3,1:NT-1)); hold on
plot(est_meas2(3,1:NT-1))
legend(["meas","est"])
title("b0z")


subplot(9,1,4)
plot(Y_noisy(4,1:NT-1)); hold on
plot(est_meas2(4,1:NT-1))
legend(["meas","est"])
title("b1x")

subplot(9,1,5)
plot(Y_noisy(5,1:NT-1)); hold on
plot(est_meas2(5,1:NT-1))
legend(["meas","est"])
title("b1y")

subplot(9,1,6)
plot(Y_noisy(6,1:NT-1)); hold on
plot(est_meas2(6,1:NT-1))
legend(["meas","est"])
title("b1z")


subplot(9,1,7)
plot(Y_noisy(7,1:NT-1)); hold on
plot(est_meas2(7,1:NT-1))
legend(["meas","est"])
title("b2x")

subplot(9,1,8)
plot(Y_noisy(8,1:NT-1)); hold on
plot(est_meas2(8,1:NT-1))
legend(["meas","est"])
title("b2y")

subplot(9,1,9)
plot(Y_noisy(9,1:NT-1)); hold on
plot(est_meas2(9,1:NT-1))
legend(["meas","est"])
title("b2z")
%%
state_est=state_est+xlp;
figure(101)
plot(state_est(8,:))



%%
function data = getParams()
    persistent params

    if isempty(params)
        parameters2;
    end
    data = params;
end