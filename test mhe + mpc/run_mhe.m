clc, clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")
import casadi.*

%load("mpc_workspace.mat")
run("run_mpc.m")
% MHE SETUP
run("mhe_setup.m")


%       MHE Simulation
%___________________________

%First we initialize MHE with initial guesses for x0 and u0
U0 = zeros(N_MHE,2);
X0 = zeros(N_MHE+1, 3);
% To improve initial guess, we use the calculated control input for U0 and
% the measured x,y for X0. Leave the initial guess for omega to be zero, no
% measurement of this
U0 = u_cl(1:N_MHE,:);
X0(:,1:2) = [y_measurements(1:N_MHE+1,1).* cos(y_measurements(1:N_MHE+1,2)), ... % [r*cos(alpha), r*sin(alpha)]
    y_measurements(1:N_MHE+1,1).* sin(y_measurements(1:N_MHE+1,2))];
mheIter=0;
st_est=[];
u_est=[];
tic

for k= 1:size(y_measurements,1) - N_MHE
mheIter=mheIter+1;

% Get the measurement window and put it as parameters in p.
args.p = [y_measurements(k:k+N_MHE,:)', u_cl(k:k+N_MHE-1,:)'];
% pass initial values
args.x0 = [reshape(X0', 3*(N_MHE+1),1); reshape(U0',2*N_MHE,1)];
sol = solver ('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
    'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
%Extract solution (estimates of applied control input and states)
U_sol = reshape(full(sol.x(3*(N_MHE+1)+1:end))',2,N_MHE)';
X_sol = reshape(full(sol.x(1:3*(N_MHE+1)))', 3, N_MHE+1)';
st_est=[st_est;X_sol(N_MHE+1,:)];
u_est =[u_est; U_sol(N_MHE,:)];
end
mhe_iter_time=toc/mheIter

save("mhe_workspace")


plotting(N_MHE, mheIter, t,t_stop, xx, u_cl, st_est, u_est, con)
Draw_MPC_PS_Obstacles (t,xx,xx1,u_cl,xs,N,robDiam,obsX, obsY, obsDiam)

