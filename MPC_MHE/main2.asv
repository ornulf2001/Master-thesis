clc,clear

addpath (genpath("C:\Users\mikae\OneDrive\Dokumenter\Industriell kyb høst 2024\Master Forprosjekt\MAGLEVREAL\O7\O7_simulering\Utlevering\simulator"))
addpath('C:\Users\mikae\OneDrive\Dokumenter\Industriell Kyb vår 2025\Master\Ny mappe\linMPCQP\ClassMPC\MPCCLASS_2DREAL\NewModel_2D')

%Dynamics
params=getParams();
% z_eq=0.0365;
% Xeq = [0; z_eq; 0; 0; 0; 0];

index = @(A, i) A(i);
fz = @(z) index(f([0, z, zeros(1, 4)]', [0, 0]', params), 5);
zeq = fzero(fz, 0.1);

Xeq = [0, zeq, zeros(1, 4)]';
Ueq = [0, 0]';
% 
% f = @(x,u) maglevSystemDynamics2d(x,u,params);
% h = @(x,u) maglevSystemMeasurements2d(x,u,params);

% [Ac, Bc, ~, ~] = finiteDifferenceLinearization(f, h, Xeq, Ueq, delta);

[Ac, Bc, Cc] = linearizeModel(@f, @h, Xeq, Ueq, params);

C = [eye(3),zeros(3,3)];
D=0;
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Tuning
X0=[0;0;0;0;0;0];
NT=200;

N_MHE=18;
N_MPC=40;
dt=0.002;

R_MHE=1*dt*diag([0.1,0.1,0.1]);
Q_MHE=2*diag([1,1,1,1,1,1]);
M_MHE=1*diag([0.1,1,0.1,2,2,2]);
noise_cov=0.05;

Q_MPC = diag([50 100 50 0.1 100 0.3]);
QN_MPC =300*Q_MPC;
R_MPC = diag([0.01 0.01]);

%Bounds
run("mpc_bounds.m")

%Run
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','off', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, QN_MPC, nStates, nControls,MPC_options);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT);
MHE_est = zeros(nStates, NT);
yCurrent=zeros(nMeasurements,NT+1);  
yCurrent(:,1)= C*X0;
newY=yCurrent(:,1);
newU=[0;0];
MPC_Xopt(:, 1) = X0;
xCurrent = X0;

for k=1:NT
    k
    
    mhe=mhe.runMHE(newY,newU);
    xEst=mhe.xCurrent;
    MHE_est(:,k)=mhe.xCurrent;
    if k<20
        U=[0;0];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xCurrent = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    elseif k>20 && k<30
        U=[1;10];
        MPC_Xopt(:,k+1) = mpc.A*MPC_Xopt(:,k) + mpc.B*U;
        xCurrent = MPC_Xopt(:,k+1);
        MPC_Uopt(:,k) = U;
        newU=MPC_Uopt(:,k);

    else
        %xCurrent-xEst
        [xNext, Uopt]=mpc.runMPC(xEst);
        MPC_Xopt(:,k+1) = xNext;
        MPC_Uopt(:,k)= Uopt;
        xCurrent=xNext;
        newU=MPC_Uopt(:,k);
    end

    yCurrent(:,k+1) = C*xCurrent + noise_cov*dt*rand(size(C*xCurrent));
    newY=yCurrent(:,k+1);
end


%Plot
figure(4)
plot(MHE_est(1,1:end)); hold on
plot(MPC_Xopt(1, 1:end));
legend('est', 'xopt')

figure(5)
plot(MHE_est(2,1:end)); hold on
plot(MPC_Xopt(2, 1:end));
legend('est', 'xopt')

figure(6)
plot(MHE_est(3,1:end)); hold on
plot(MPC_Xopt(3, 1:end));
legend('est', 'xopt')
% mpc.plotResults(MPC_Xopt,MPC_Uopt)



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
A = round(A, 5);

for i = 1:nInputs
    B(:,i) = (f(x_eq,u_eq+delta*(i==1:nInputs)')-f(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
B = round(B, 5);

for i = 1:nStates
    C(:,i) = (h(x_eq+delta*(i==1:nStates)',u_eq)-h(x_eq-delta*(i==1:nStates)',u_eq))/(2*delta);
end
C = round(C,5);

for i = 1:nInputs
    D(:,i) = (h(x_eq,u_eq+delta*(i==1:nInputs)')-h(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
D = round(D,5);
end

function [A,B,C] = linearizeModel(f,h,xlp,ulp,params)
    delta = 1e-10;
    
    A = zeros(6, 6);    
    for i = 1:6
        xp = xlp;
        xn = xlp;
    
        xp(i) = xp(i) + delta;
        xn(i) = xn(i) - delta;
        
        dxp = f(xp, ulp, params);
        dxn = f(xn, ulp, params);
        
        A(:, i) = (dxp - dxn)/(2*delta);
    end
    
    B = zeros(6, 2);    
    for i = 1:2
        up = ulp;
        un = ulp;
    
        up(i) = up(i) + delta;
        un(i) = un(i) - delta;
        
        dxp = f(xlp, up, params);
        dxn = f(xlp, un, params);
        
        B(:, i) = (dxp - dxn)/(2*delta);
    end
    
    C = zeros(2, 6);    
    for i = 1:6
        xp = xlp;
        xn = xlp;
    
        xp(i) = xp(i) + delta;
        xn(i) = xn(i) - delta;
        
        yp = h(xp, params);
        yn = h(xn, params);
        
        C(:, i) = (yp - yn)/(2*delta);
    end
end
