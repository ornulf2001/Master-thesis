clc;
clear all;
run('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\system_parameters\parameters.m');
addpath('C:\Program Files\casadi-3.6.7-windows64-matlab2018b');
run('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\system_parameters\parameters.m')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\model_implementations\fields_and_force')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\model_implementations\dynamics_and_measurements')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\utilities\general')

c = @(x,u) maglevSystemDynamics2d(x,u,params);
h = @(x,u) maglevSystemMeasurements2d(x,u,params);

import casadi.*

Rad = 0.04;
M = 0.05;
J = 1.1;
m = 1;
gr = 9.81;
mu0 = 4*pi*1e-7;

z_eq = 0.365;
delta = 1e-5;              % Perturbasjon fra lineariseringspunktet
x_eq = [0,z_eq,0,0,0,0]';
u_eq = [0,0]';

x = SX.sym('x'); z = SX.sym('z'); theta = SX.sym('theta');
xdot = SX.sym('xdot'); zdot = SX.sym('zdot'); thetadot = SX.sym('thetadot');
states = [x;z;theta; xdot; zdot; thetadot]; n_states = length(states);

ux = SX.sym('ux'); uz = SX.sym('uz');
controls = [ux;uz]; n_controls = length(controls);

rn = [x - Rad*cos(theta);
      0;
      z - Rad*sin(theta)];

rp = [x + Rad*cos(theta);
      0;
      z + Rad*sin(theta)];

Ip = (ux + uz)/2;
In = (uz - ux)/2;


mn = (m + mu0*In).*[sin(theta); 0; -cos(theta)]';
mp = (m + mu0*Ip).*[sin(theta); 0; -cos(theta)]';
ml = [0; 0; m];

norm_rp = (rp(1)^2 + rp(3)^2)^(1/2);
norm_rn = (rn(1)^2 + rn(3)^2)^(1/2);

phip = mu0/(4*pi) * ((mp*rp/(norm_rp^3)));
phin = mu0/(4*pi) * ((mn*rn)/norm_rn^3);
phi = phip + phin;

dphi_dx = jacobian(phi, x);
d2phi_dxdz = jacobian(dphi_dx, z);
dphi_dz = jacobian(phi, z);
d2phi_dz2 = jacobian(dphi_dz, z);

xddot = (1/M) * (-m*d2phi_dxdz);
zddot = (-(m/M*d2phi_dz2));
thetaddot = (1/J) * (-m*dphi_dx);

h = 0.1; % sampling time [s]
N = 3; % prediction horizon
% rob_diam = 0.3;

ux_max = 1; ux_min = -ux_max;
uz_max = 1; uz_min = -uz_max;

rhs = [xdot; zdot; thetadot; xddot; zddot; thetaddot]; % system r.h.s


f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X',n_states,(N+1));
% A Matrix that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = diag([10000 10000 100 10 10 10]);
R = 0.1*diag([0.1 0.1]);

st = X(:, 1);
g = [g; st - P(1:6)];


for k = 1:N
    st = X(:, k);
    con = U(:, k);
    obj = obj + (st - P(7:12))' * Q * (st - P(7:12)) + con'*R*con;
    st_next = X(:, k+1);
    k1 = f(st, con);
    k2 = f(st + h/2*k1, con);
    k3 = f(st + h/2*k2, con);
    k4 = f(st + h*k3, con);
    stNextRK4 = st + h/6*(k1 + 2*k2 + 2*k3 + k4);
    g = [g; st_next - stNextRK4];
end



% make the decision variables one column vector
OPT_variables = [reshape(X, 6*(N+1), 1); reshape(U, 2*N, 1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;

args.lbg(1:6*(N+1)) = 0;
args.ubg(1:6*(N+1)) = 0;

args.lbx(1:6:6*(N+1), 1) = 0;
args.ubx(1:6:6*(N+1), 1) = 1;
args.lbx(2:6:6*(N+1), 1) = 0;
args.ubx(2:6:6*(N+1), 1) = 10;
args.lbx(3:6:6*(N+1), 1) = -pi/4;
args.ubx(3:6:6*(N+1), 1) = pi/4;

args.lbx(6*(N+1) + 1:2:6*(N+1)+2*N, 1) = ux_min;
args.ubx(6*(N+1) + 1:2:6*(N+1)+2*N, 1) = ux_max;
args.lbx(6*(N+1) + 2:2:6*(N+1)+2*N, 1) = uz_min;
args.ubx(6*(N+1) + 2:2:6*(N+1)+2*N, 1) = uz_max;
% inequality constraints (state constraints)
% args.lbg = repmat([-1; 0; -pi/2], N+1, 1);  % lower bound of the states x and y
% args.ubg = repmat([1; 1; pi/2], N+1, 1);   % upper bound of the states x and y 
% 
% 
% % input constraints
% args.lbx(1:2:2*N-1,1) = ux_min; args.lbx(2:2:2*N,1)   = uz_min;
% args.ubx(1:2:2*N-1,1) = ux_max; args.ubx(2:2:2*N,1)   = uz_max;


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0; 0.3; 0; 0; 0; 0];    % initial condition.
xs = [0.01; 0.2; 0; 0; 0; 0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

X0 = repmat(x0, 1, N+1)';

u0 = zeros(N,2);  % two control inputs 

sim_tim = 30; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];


% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-3 && mpciter < sim_tim / h)

    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = [reshape(X0', 6*(N+1), 1); reshape(u0', 2*N, 1)];
    
    % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x(6*(N+1)+1:end))', 2, N)';
    % ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:, 1:6, mpciter+1) = reshape(full(sol.x(1:6*(N+1)))', 6, N+1)';

    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = shift(h, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0;  
    X0 = reshape(full(sol.x(1:6*(N+1)))', 6, N+1)';
    X0 = [X0(2:end, :); X0(end, :)];
    mpciter
    mpciter = mpciter + 1;
end
toc

% Plotting x1 and x2
t = 0:mpciter;

figure(1);
subplot(3, 1, 1);
plot(t, xx(1, :), 'k',  'LineWidth', 0.5);
grid on;
hold on;
legend('$x_1$', 'Interpreter', 'latex')
title('States x_1 and x_2');
ylabel('$x_1$', 'Interpreter', 'latex')
subplot(3, 1, 2);
plot(t, xx(2, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_2$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
subplot(3, 1, 3);
plot(t, xx(3, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_3$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_3$', 'Interpreter', 'latex')

figure(2);
subplot(3, 1, 1);
plot(t, xx(4, :), 'k',  'LineWidth', 0.5);
grid on;
hold on;
legend('$x_1$', 'Interpreter', 'latex')
title('States x_1 and x_2');
ylabel('$x_1$', 'Interpreter', 'latex')
subplot(3, 1, 2);
plot(t, xx(5, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_2$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
subplot(3, 1, 3);
plot(t, xx(6, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_3$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_3$', 'Interpreter', 'latex')

% Plotting u as stairs
figure(3);
subplot(2, 1, 1);
stairs(t(1:end-1), u_cl(:,1), 'r.-',  'LineWidth', 0.5);
grid on;
hold on;
subplot(2, 1, 2)
stairs(t(1:end-1), u_cl(:,2), 'r.-',  'LineWidth', 0.5);
grid on;
legend('$u$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
xlabel('$Time [s]$', 'Interpreter', 'latex');

function [t0, x0, u0] = shift(T, t0, x0, u, f)
st = x0;
con = u(1, :)';
f_val = f(st, con);
st = st + (T*f_val);
x0 = full(st);

t0 = t0 + T;
u0 = [u(2:size(u,1), :); u(size(u, 1), :)];
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
