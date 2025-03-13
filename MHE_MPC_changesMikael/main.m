clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*

params=getParams();
z_eq=0.0365;
Xeq = [0; z_eq; 0; 0; 0; 0];
Ueq=[0;0];
delta = 1e-6;
f = @(x,u) maglevSystemDynamics2d(x,u,params);
h = @(x,u) maglevSystemMeasurements2d(x,u,params);
[Ac, Bc, ~, ~] = finiteDifferenceLinearization(f, h, Xeq, Ueq, delta);
C = [eye(3),zeros(3,3)];
D=0;

%Run msd simulation
dt=0.002;
X0=[0.001;0.003;00.003;0;0;0];
%[z0block,C,D,T_sim,U_list,X_sim,Y_sim,Y_noisy]=sim_2d(dt,X0,T_END,params,Ac,Bc);
%%%%%%%% Setup %%%%%%%%%

% Sizes
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);
NT=500;

% Weights and 
N_MHE=12;
N_MPC=80;
R_MHE=2*dt*diag([0.3,0.8,1]);
Q_MHE=diag([1,1,1,1,1,1]);
M_MHE=diag([1,1,1,3,3,3]);

Q_MPC = diag([20 10 20 0.5 0.5 0.5]);
QN_MPC =30*Q_MPC;
R_MPC = diag([0.001 0.001]);

%Bounds MPC
lbx = -inf;
ubx = inf;
lbz = -inf;
ubz = inf;
lbtheta = -inf;
ubtheta = inf;

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
lbU = repmat(lbU1, N_MPC, 1);
ubU = repmat(ubU1, N_MPC, 1);

lb = [lbX; lbU];
ub = [ubX; ubU];



%%%%%%%%% Run  %%%%%%%%%%
MHE_options = optimoptions("quadprog","Display","off", "Algorithm","interior-point-convex");
mhe = MHEclass(N_MHE,Ac,Bc,C,Q_MHE,R_MHE,M_MHE,X0,dt,MHE_options);

MPC_options = optimset('Display','on', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
mpc = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q_MPC, R_MPC, QN_MPC, nStates, nControls,MPC_options);
MPC_Xopt = zeros(nStates, NT);
MPC_Uopt = zeros(nControls, NT);
MPC_Xopt(:, 1) = X0;

xCurrent=X0;
x_ests=[];
yCurrent=zeros(size(C*X0,1),NT+1);   
yCurrent(:,1) = C*X0;

Y_noisy=[]
U_list=[]
    % Run class version of MHE
    
    for k=2:NT
    
        
        %if ~mhe.isReadyToRun
            if k < 2*mhe.yBufferCount+2
            [Xopt, Uopt]=mpc.runMPC(xCurrent);
            xEst=xCurrent;
       % else
            else 
             %   xEst=mhe.A*MPC_Xopt(:,k)+mhe.B*MPC_Uopt(:,k);
              [Xopt, Uopt]=mpc.runMPC(xEst);
            end
        
        MPC_Uopt(:,k-1)= Uopt;
        MPC_Xopt(:,k) = Xopt;

        xCurrent=Xopt;
        x_sim_next=mhe.A*MPC_Xopt(:,k-1)+mhe.B*Uopt;   %maglevSystemDynamics2d(MPC_Xopt(:,k-1),Uopt,params);
        yCurrent(:,k) = C*x_sim_next + 0*dt*rand(size(C*x_sim_next));
        newY=yCurrent(:,k);
        newU=Uopt;
        if ~mhe.isReadyToRun
            mhe=mhe.bufferInitialData(newY, newU);
            
        else
            mhe=mhe.runMHE(newY,newU);
            xEst=mhe.xCurrent;
            x_ests=[x_ests,xEst];
        end
        
        

    end   


    % mhe=mhe.reset(X0);
    % while ~mhe.isReadyToRun
    %    newY=yCurrent(:,mhe.yBufferCount);
    %    newU=MPC_Uopt(:,mhe.uBufferCount);
    %    mhe=mhe.bufferInitialData(newY, newU);
    % end
    % xsol2=[];
    % xsol2(:,1)=X0;
    % for k=2:size(yCurrent,2)-1
    %     newY=yCurrent(:,k);
    %     newU=MPC_Uopt(:,k-1);
    %     mhe=mhe.runMHE(newY,newU);
    %     mhe.xCurrent
    %     xsol2(:,k)=mhe.xCurrent;
    % end
    % plot(xsol2(2,:));hold on
    % plot(MPC_Xopt(2,N_MHE+1:end)); 
    %plot(yCurrent(2,1:end-1))
    figure(1)
    plot(MPC_Xopt(1,N_MHE+1:end)); hold on
    plot(x_ests(1,:))

    figure(2)
    plot(MPC_Xopt(2,N_MHE+1:end)); hold on
    plot(x_ests(2,:))

    figure(3)
    plot(MPC_Xopt(3,N_MHE+1:end)); hold on
    plot(x_ests(3,:))


    %ylim([-0.1,0.1])
    %plot(yCurrent(2,1:end))
    %legend(["mpc opt"]);%,"meas"])
    %plot([zeros(1,mhe.yBufferCount),x_ests(2,:)])

     


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
    zBlock(5) = Az1*z_eq;

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
A = round(A,5);

for i = 1:nInputs
    B(:,i) = (f(x_eq,u_eq+delta*(i==1:nInputs)')-f(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
B = round(B,5);

for i = 1:nStates
    C(:,i) = (h(x_eq+delta*(i==1:nStates)',u_eq)-h(x_eq-delta*(i==1:nStates)',u_eq))/(2*delta);
end
C = round(C,5);

for i = 1:nInputs
    D(:,i) = (h(x_eq,u_eq+delta*(i==1:nInputs)')-h(x_eq,u_eq-delta*(i==1:nInputs)'))/(2*delta);
end
D = round(D,5);
end