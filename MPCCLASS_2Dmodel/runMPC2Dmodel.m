clc, clear all
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell Kyb vår 2025\Master\Ny mappe\qpOASES-master\qpOASES-master\interfaces\matlab')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\system_parameters')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\model_implementations\fields_and_force')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\model_implementations\dynamics_and_measurements')
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator\utilities\general')

params=getParams()


f = @(x,u) maglevSystemDynamics2d(x,u,params);
h = @(x,u) maglevSystemMeasurements2d(x,u,params);

% RunMPC for 2D-maglev model

% params
mu0 = params.physical.mu0;
gr = params.physical.g;
Rad = params.solenoids.x(1); %Avstand fra nullpunkt til solenoid
rl = params.magnet.r; %Radius svevemagnet
rp = params.permanent.r(1); %Radius permanent magnet
ll = params.magnet.l; %Lengde svevemagnet
lp = params.permanent.l(1); %Lengde permanent magnet
M = params.magnet.m; %Mass
J = params.magnet.I(1); %Moment of inertia
ml = abs(params.magnet.J)*pi*rl^2*ll/mu0; % magnetisk moment svevemagnet
m = abs(params.permanent.J)*pi*rp^2*lp/mu0; % magnetisk moment permanent magnet
% 

z_eq = 0.0365;
dt = 0.01;
X0 = [-0.05; 0.2; pi/8; 0.01; 0.05; -0.1];
xRef = [0; z_eq; pi/12; 0; 0; 0];
% zBlock = zeros(size(X0));
N_MPC = 2;
NT = 120;

Q = diag([3 1 1 2 2 2]);
QN = Q;
R = diag([1 1]);

Xeq = [0; z_eq; 0; 0; 0; 0];
UeqGuess = [0; 0];
delta = 1e-5;

nStates = size(Q, 2);
nControls = size(R, 2);


%Bounds on states and inputs
lbx = -0.1;
ubx = 0.1;
lbz = -1;
ubz = 1;
lbtheta = -pi/4;
ubtheta = pi/4;

lbxdot = -inf;
ubxdot = inf;
lbzdot = -inf;
ubzdot = inf;
lbthetadot = -inf;
ubthetadot = inf;

lbux = -inf;
ubux = inf;
lbuz = -inf;
ubuz = inf;

lbX1 = [lbx; lbz; lbtheta; lbxdot; lbzdot; lbthetadot];
ubX1 = [ubx; ubz; ubtheta; ubxdot; ubzdot; ubthetadot];
lbX = repmat(lbX1, N_MPC+1, 1);
ubX = repmat(ubX1, N_MPC+1, 1);

lbU1 = [lbux; lbuz];
ubU1 = [ubux; ubuz];
lbU = repmat(lbU1, N_MPC+1, 1);
ubU = repmat(ubU1, N_MPC+1, 1);

lb = [lbX; lbU];
ub = [ubX; ubU];

[Ac, Bc, zBlock] = linMat(z_eq, params, nStates);

option = optimset('Display', 'iter', 'TolFun', 1e-6);
Ueq = fsolve(@(u) maglevSystemDynamics2d(Xeq, u, params), UeqGuess, option);

% [Ac, Bc, ~, ~] = finiteDifferenceLinearization(f, h, Xeq, Ueq, delta);


mpc = MPCclass(N_MPC, Ac, Bc, xRef, zBlock, X0, dt, lb, ub, Q, R, QN, NT, nStates, nControls);

[Xhis, Uopt] = mpc.runMPC();
mpc.plotResults(Xhis, Uopt);










%% System dynamics for 2D system
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



function [A, B, zBlock] = linMat(z_eq, params, nStates)
    

    %params
    mu0 = params.physical.mu0;
    gr = params.physical.g;
    Rad = params.solenoids.x(1); %Avstand fra nullpunkt til solenoid
    rl = params.magnet.r; %Radius svevemagnet
    rp = params.permanent.r(1); %Radius permanent magnet
    ll = params.magnet.l; %Lengde svevemagnet
    lp = params.permanent.l(1); %Lengde permanent magnet
    M = params.magnet.m; %Mass
    J = params.magnet.I(1); %Moment of inertia
    ml = abs(params.magnet.J)*pi*rl^2*ll/mu0; % magnetisk moment svevemagnet
    m = abs(params.permanent.J)*pi*rp^2*lp/mu0; % magnetisk moment permanent magnet

    nw = params.solenoids.nw;
    rs = params.solenoids.nw;
    As = pi*rs^2;

    % Koeffisienter for lineær modell
    q  = Rad^2 + z_eq^2;
    
    p1 = -(4*Rad^4 - 27*Rad^2*z_eq^2 + 4*z_eq^4)/2;
    p2 = -(z_eq*(11*Rad^4 - 23*Rad^2*z_eq^2 + z_eq^4))/4;
    
    p3 = -10*(Rad*(Rad^2 - 4*z_eq^2))/2;
    
    p4 = 2*(-3*Rad^4 + 24*Rad^2*z_eq^2 - 8*z_eq^4);
    
    p5 = 1e3*(z_eq*(3*Rad^2 - 2*z_eq^2))/2;
    
    p6 = 20*(z_eq*(4*Rad^2 - z_eq^2));
    p7 = 20*(-Rad^4 + 13*Rad^2*z_eq^2 - z_eq^4);
    
    p8 = 200*(Rad*z_eq)/2;
    % 
    % Ax1 = (3*mu0*m^2)/(2*pi*M*q^(9/2))*p1;
    % Ax2 = (3*mu0*m^2)/(2*pi*M*q^(9/2))*p2;
    % Bx = -(3*mu0^2*m)/(2*pi*M*q^(7/2))*p3;
    % 
    % Az1 = -(3*mu0*m^2)/(2*pi*M*q^(9/2))*p4;
    % Bz = -(3*mu0^2*m)/(2*pi*M*q^(7/2))*p5;
    % 
    % Ath1 = (3*mu0*m^2)/(2*pi*J*q^(7/2))*p6;
    % Ath2 = -(3*mu0*m^2)/(2*pi*J*q^(7/2))*p7;
    % Bth = (3*mu0^2*m)/(2*pi*J*q^(5/2))*p8;

    Ax1 = (3*mu0*ml*m*p1)/(2*pi*M*q^(9/2));
    Ax2 = (3*mu0*ml*m*p2)/(2*pi*M*q^(9/2));
    Bx = -(3*mu0*ml*m*nw*As*p3)/(2*pi*M*m*q^(7/2));

    Az1 = (-3*mu0*ml*m*p4)/(2*pi*M*q^(9/2));
    Bz = (-3*mu0*ml*m*nw*As*p5)/(2*pi*M*m*q^(7/2));

    Ath1 = (3*mu0*ml*m*p6)/(2*pi*J*q^(7/2));
    Ath2 = (-3*mu0*ml*m*p7)/(2*pi*J*q^(7/2));
    Bth = (3*mu0*ml*m*nw*As*p8)/(2*pi*J*m*q^(5/2));
    
    A = [0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        Ax1 0 Ax2 0 0 0;
        0 Az1 0 0 0 0;
        Ath1 0 Ath2 0 0 0];
    
    B = [0 0;
        0 0;
        0 0;
        Bx 0;
        0 Bz;
        Bth 0];

    zBlock = zeros(nStates, 1);
    zBlock(5) = Az1*-z_eq;

end



function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end

function [A,B,C,D] = finiteDifferenceLinearization(f,h,x_eq,u_eq,delta)
nStates     = length(x_eq);
nInputs     = length(u_eq);
nOutputs    = length(h(x_eq,u_eq));

A = zeros(nStates);
B = zeros(nStates,nInputs);
C = zeros(nOutputs,nStates);
D = zeros(nOutputs,nInputs);

for i = 1:nStates
    A(:,i) = (f(x_eq+delta*(i==1:nStates)',u_eq)-f(x_eq-delta*(i==1:nStates)',u_eq))/(2*delta);
end
A = A

for i = 1:nInputs
    B(:,i) = (f(x_eq,u_eq+delta*(i==1:nInputs)')-f(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
B = B

for i = 1:nStates
    C(:,i) = (h(x_eq+delta*(i==1:nStates)',u_eq)-h(x_eq-delta*(i==1:nStates)',u_eq))/(2*delta);
end
C = round(C,5);

for i = 1:nInputs
    D(:,i) = (h(x_eq,u_eq+delta*(i==1:nInputs)')-h(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
D = round(D,5);
end