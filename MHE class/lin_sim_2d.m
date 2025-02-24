clc,clear

dt=0.005;
T_sim = 0:dt:5; 
N_MHE=5;

z0=0.0365;
params=getParams();
[Ac,Bc,~]=f_func_lin(z0,params);

z0block=zeros(size(Ac,1),1);

A=expm(Ac*dt);
B=inv(Ac) * (A - eye(size(Ac))) * Bc;
C = [eye(2),zeros(2,4)];

x0_sim=[0.01;0.06;0.001;0;0;0];
u0 = [0;0];
u=u0;
%u_sim=u0;
xref=[0;z0];
%X_sim = zeros(length(x0_sim), length(T_sim)); % State storage
X_sim = x0_sim; % Initial state
X_est=x0_sim;
%u = zeros(size(B,2),length(T_sim));

Qk = diag([10,10,10,1,10,1]);
Rk = diag([0.01,0.01]);
[K,~,~]  = dlqr(A,B,Qk,Rk);
F = inv(C * ((Bc*K - Ac) \ Bc));

QL = diag([10,10,10,1,10,1]);
RL = diag([0.01,0.01]);
%L_obs=dlqr(A',C',QL,RL)';

poles=[ 0.8, 0.85, 0.82, 0.83, 0.89, 0.86];
L_obs=place(A',C',poles)';

function yk=measure(xk,C)
    yk = C*xk + 0*randn(size(C*xk,1),1);
end


R=diag([0.001, 0.001, 0.001]);
Q=diag([1,1,1,1,1,1]);
M=diag([1,1,1,1,1,1]);
mhe=MHEclass_linear(N_MHE,Ac,Bc,C,Q,R,M,z0block,x0_sim,dt);

for i=1:100
    X_sim = [X_sim, A*X_sim(:,end) + B*u(:,end)];
        if ~mhe.isReadyToRun
            newY = measure(X_sim(:,end), C);
            correction = L_obs*(newY - C*X_est(:,end));

            new_est = A*X_est(:,end) + B*u(:,end) + 0*correction
            X_est = [X_est, new_est];

            newU = -K*X_est(:,end) + F*xref;  
            %newU_sim=-K*(X_sim(:,end))+ F*xref;
            u=[u,newU];
            %u_sim=[u_sim,newU_sim];
            %mhe=mhe.bufferInitialData(newY, newU);
        end
        %newY=measure(X_sim(:,end),C);


end
PPP=mhe.P;
% for k=1:length(T_sim)
%     uk = -K*(X_sim(:,k))+ F*xref;
%     X_sim(:, k+1)=A*X_sim(:, k)+B*uk ;
% end
% 
% plot(X_sim(2,:))
% ylim([-0.2,0.2])


function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end
function dx = maglevSystemDynamics2d(x,u,params)
    % ### HACK from 2d to 3d ###
    x_full = [
        x(1),0,x(2),0,x(3),0,...
        x(4),0,x(5),0,x(6),0
        ]';
    u_full = [u(1), u(2), 0, 0]';
    % ##########################
    
    dx = maglevSystemDynamics(x_full,u_full,params);
    
    % ### HACK from 3d to 2d ###
    dx = dx([1,3,5,7,9,11]');
    % ##########################
end

function y = maglevSystemMeasurements2d(x,u,params)
    % ### HACK from 2d to 3d ###
    x_full = [
        x(1),0,x(2),0,x(3),0,...
        x(4),0,x(5),0,x(6),0
        ]';
    u_full = [u(1), u(2), 0, 0]';
    % ##########################
    
    y = maglevSystemMeasurements(x_full,u_full,params);
    
    % ### HACK from 3d to 2d ###
    y = y([1,3]');
    % ##########################
end