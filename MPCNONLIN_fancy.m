


clc;
clear all;

%run('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\system_parameters\parameters.m')
addpath('C:\Users\ornul\AppData\Roaming\MathWorks\CasADi');


% Importing CasADi toolbox
import casadi.*


% Chosing parameter values
Rad = 0.025;
M = 0.06;
J = 1;
m = 2.47;
gr = 9.81;
mu0 = 4*pi*1e-7;

% Rad = params.r;
% M = params.m;
% J = params.J;
% m = params.
% gr = params.g;
% mu0 = params.mu0;

% SX symbolic expression for states
x = SX.sym('x'); 
z = SX.sym('z'); 
theta = SX.sym('theta');
xdot = SX.sym('xdot'); 
zdot = SX.sym('zdot'); 
thetadot = SX.sym('thetadot');
states = [x; z; theta; xdot; zdot; thetadot]; 
n_states = length(states);

% SX symbolic expression for inputs
ux = SX.sym('ux'); 
uz = SX.sym('uz');
controls = [ux;uz]; 
n_controls = length(controls);


params=getParams()


% Calculating the Nonlinear dynamics
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

% System model equations
xddot = (1/M) * (-m*d2phi_dxdz);
zddot = (-(m/M*d2phi_dz2) - gr);
thetaddot = (1/J) * (-m*dphi_dx);


h = 0.2; % Sampling time
N = 3; % Prediction Horizon

% Control input bounds
ux_max = 1; 
ux_min = -ux_max;
uz_max = 1; 
uz_min = -uz_max;

% Symbolic Right Hand Side vector
rhs = [xdot; zdot; thetadot; xddot; zddot; thetaddot]; 

f = Function('f',{states,controls},{rhs}); % Nonlin mapping function F(x, u)
U = SX.sym('U',n_controls,N); % Matrix storage for control inputs
P = SX.sym('P',n_states + n_states); % Matrix storage for current states and reference state
X = SX.sym('X',n_states,(N+1)); % Matrix storage for states

obj = 0; % Objective function
g = [];  % Constraints vector

% Weight matrices Q and R
Q = diag([10 10 10 10 10 10]);
R = diag([10 10]);


st = X(:, 1); % Init state vector
g = [g; st - P(1:6)]; % Init constraint vector

% For loop to construct constraints and objective function
for k = 1:N

    st = X(:, k);
    con = U(:, k);
    obj = obj + (st - P(7:12))' * Q * (st - P(7:12)) + ...
        con' * R * con; % Objective function 

    stNext = X(:, k+1);

    k1 = f(st, con);
    k2 = f(st + h/2*k1, con);
    k3 = f(st + h/2*k2, con);
    k4 = f(st + h*k3, con);

    stRK4 = st + h/6*(k1 + 2*k2 + 2*k3 + k4); % RK4 method for calculating xi

    g = [g; stNext - stRK4]; % Constraint vector

end


% NLP solver setup
OPT_variables = [reshape(X, 6*(N+1), 1); reshape(U, 2*N, 1)]; % Reshape decition variables as vector

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P); %

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts); % Ipopt as solver


args = struct;


% Inequality constraints, bounds
args.lbg(1:6*(N+1)) = 0;
args.ubg(1:6*(N+1)) = 0;

args.lbx(1:6:6*(N+1), 1) = 0;
args.ubx(1:6:6*(N+1), 1) = 0.5;
args.lbx(2:6:6*(N+1), 1) = 0;
args.ubx(2:6:6*(N+1), 1) = 0.05;
args.lbx(3:6:6*(N+1), 1) = -pi/4;
args.ubx(3:6:6*(N+1), 1) = pi/4;
args.lbx(4:6:6*(N+1), 1) = 0;
args.ubx(4:6:6*(N+1), 1) = 1;
args.lbx(5:6:6*(N+1), 1) = 0;
args.ubx(5:6:6*(N+1), 1) = 10;
args.lbx(6:6:6*(N+1), 1) = -pi/4;
args.ubx(6:6:6*(N+1), 1) = pi/4;

args.lbx(6*(N+1) + 1:2:6*(N+1)+2*N, 1) = ux_min;
args.ubx(6*(N+1) + 1:2:6*(N+1)+2*N, 1) = ux_max;
args.lbx(6*(N+1) + 2:2:6*(N+1)+2*N, 1) = uz_min;
args.ubx(6*(N+1) + 2:2:6*(N+1)+2*N, 1) = uz_max;


% MPC simulation setup

t0 = 0;
x0 = [0; 0.03; 0; 0; 0; 0];    % initial point
xRef = [0; 0.04; 0; 0; 0; 0]; % Reference point

Xhis(:, 1) = x0; % Preallocating Xhis to store history of states
t(1) = t0;

X0 = repmat(x0, 1, N+1)';

u0 = zeros(N, 2);  % Initial control inputs as zero

simTime = 10; % Simulation time


it = 0;
uopt = [];
init_plot(x0)


% If the error between x0 and xRef is lesser than 1e-3 it doesnt work
%NMPCloop = tic;
while (norm((x0 - xRef), 2) > 1e-3 && it < simTime / h)

    args.p   = [x0; xRef]; % Parameter vector
    args.x0 = [reshape(X0', 6*(N+1), 1); reshape(u0', 2*N, 1)]; 
    
    
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p); % NLP solver Ipopt
    
    u = reshape(full(sol.x(6 * (N+1)+1:end))', 2, N)'; % Extract the optimal control input u
      
    uopt = [uopt ; u(1,:)]; % Store updated control input

    t(it + 1) = t0;
    
    % Using a premade function to get init for next opt step
    [t0, x0, u0] = UpdateINIT(h, t0, x0, u, maglevSystemDynamics2d(x0,u(1, :)',params));
    
    Xhis(:, it+2) = x0; % Updating Xhis as storage for each step
    X0 = reshape(full(sol.x(1:6*(N+1)))', 6, N+1)';
    X0 = [X0(2:end, :); X0(end, :)];
    %it % Printing the number of iterations
    it = it + 1;
    update_plot(x0)
end
%toc


% Plotting results

t = 0:it;

figure(1);
subplot(3, 1, 1);
plot(t, Xhis(1, :), 'k',  'LineWidth', 0.5);
grid on;
hold on;
legend('$x$', 'Interpreter', 'latex')
title('States');
ylabel('$x$', 'Interpreter', 'latex')
subplot(3, 1, 2);
plot(t, Xhis(2, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$z$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$z$', 'Interpreter', 'latex')
subplot(3, 1, 3);
plot(t, Xhis(3, :), 'g',  'LineWidth', 0.5);
grid on;
legend('$\theta$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$\theta$', 'Interpreter', 'latex')

figure(2);
subplot(3, 1, 1);
plot(t, Xhis(4, :), 'k',  'LineWidth', 0.5);
grid on;
hold on;
legend('$x_1$', 'Interpreter', 'latex')
title('States x_1 and x_2');
ylabel('$x_1$', 'Interpreter', 'latex')
subplot(3, 1, 2);
plot(t, Xhis(5, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_2$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
subplot(3, 1, 3);
plot(t, Xhis(6, :), 'b',  'LineWidth', 0.5);
grid on;
legend('$x_3$', 'Interpreter', 'latex')
xlabel('$Time [s]$', 'Interpreter', 'latex')
ylabel('$x_3$', 'Interpreter', 'latex')

% Plotting u as stairs
figure(3);
subplot(2, 1, 1);
stairs(t(1:end-1), uopt(:,1), 'r.-',  'LineWidth', 0.5);
legend('$ux$', 'Interpreter', 'latex');
ylabel('$ux$', 'Interpreter', 'latex');
title('Control Inputs')
grid on;
hold on;
subplot(2, 1, 2)
stairs(t(1:end-1), uopt(:,2), 'r.-',  'LineWidth', 0.5);
grid on;
legend('$uz$', 'Interpreter', 'latex');
ylabel('$uz$', 'Interpreter', 'latex');
xlabel('$Time [s]$', 'Interpreter', 'latex');



% Function to get the initial for next opt step
function [t0, x0, u0] = UpdateINIT(T, t0, x0, u, f)
st = x0;
%con = u(1, :)';

fval = f;
st = st + (T * fval);
x0 = full(st);

t0 = t0 + T;
u0 = [u(2:size(u, 1), :); u(size(u, 1), :)];
end



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

function init_plot(x0)
    % Persistent variables for plotting
    MaglevPlotFig = getappdata(0, 'MaglevPlotFig');

    if isempty(MaglevPlotFig) || ~ishandle(MaglevPlotFig)
        % Create figure and axes
        MaglevPlotFig = figure('Name', 'Maglev Animation', 'NumberTitle', 'off');
        MaglevPlotAxes = axes('Parent', MaglevPlotFig);
        grid(MaglevPlotAxes, 'on'); hold(MaglevPlotAxes, 'on'); box(MaglevPlotAxes, 'on');
        daspect(MaglevPlotAxes, [1, 1, 1]);
        xlim(MaglevPlotAxes, [-1, 1]);
        ylim(MaglevPlotAxes, [-0.05, 1]);
        xlabel(MaglevPlotAxes, '$x$','Interpreter','latex','FontSize',14);
        ylabel(MaglevPlotAxes, '$z$','Interpreter','latex','FontSize',14);

        % Draw ground
        yline(MaglevPlotAxes, 0, 'k', 'LineWidth', 2);

        % Draw solenoids
        w = 2*0.0092;
        h = 0.0120;
        rectangle('Parent', MaglevPlotAxes, 'Position',[0.02-w/2,0,w,h], 'EdgeColor', 'k', 'FaceColor',[0.72, 0.45, 0.2], 'LineWidth',2);
        rectangle('Parent', MaglevPlotAxes, 'Position',[-0.02-w/2,0,w,h], 'EdgeColor', 'k', 'FaceColor',[0.72, 0.45, 0.2], 'LineWidth',2);

        % Draw maglev vehicle
        MaglevHandle = create_maglev(MaglevPlotAxes);
        set(MaglevHandle, 'Matrix', makehgtform('translate', [x0(1), x0(2), 0], 'zrotate', -x0(3)));

        % Store handles in application data
        setappdata(0, 'MaglevPlotFig', MaglevPlotFig);
        setappdata(0, 'MaglevHandle', MaglevHandle);
        setappdata(0, 'LastUpdateTime', 0);
    end
end

function update_plot(x)
    % Retrieve persistent data
    MaglevHandle = getappdata(0, 'MaglevHandle');

    % Update maglev position and orientation
    set(MaglevHandle, 'Matrix', makehgtform('translate', [x(1), x(2), 0], 'zrotate', x(3)));
    drawnow limitrate nocallbacks; % Efficient drawing
end

function MaglevHandle = create_maglev(ax)
    % Derived parameters
    w_body = 2*0.0250;
    h_body = 0.0050;
    r_center = 0.001;
    linewidth = 1;

    % Maglev transform object
    MaglevHandle = hgtransform('Parent', ax);

    % Body
    rectangle('Parent', MaglevHandle, 'Position', [-w_body/2, -h_body/2, w_body, h_body], ...
        'EdgeColor', 'k', 'FaceColor', [0.5, 0.5, 0.5], 'LineWidth', linewidth);

    % Center marker
    rectangle('Parent', MaglevHandle, 'Position', [-r_center/2, -r_center/2, r_center, r_center], ...
        'Curvature', [1, 1], 'EdgeColor', 'k', 'FaceColor', 'k', 'LineWidth', linewidth);
end