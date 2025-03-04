clc,clear
CasADiPath = "C:\Users\ornul\AppData\Roaming\MathWorks\CasADi";
addpath(CasADiPath)
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*



%%%% Setup %%%%%%%%%


%Dynamical model


% Sizes
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%Misc.
N_MHE=5;

R=diag([0.001, 0.001, 0.001]);
Q=diag([1,1,1,1,1,1]);
M=diag([1,1,1,1,1,1]);


%Nonlinear dynamics and measurements
% x = SX.sym('x'); 
% z = SX.sym('z'); 
% theta = SX.sym('theta');
% xdot = SX.sym('xdot'); 
% zdot = SX.sym('zdot'); 
% thetadot = SX.sym('thetadot');
% states = [x; z; theta; xdot; zdot; thetadot]; 
% %nStates = length(states);
% 
% ux = SX.sym('ux'); 
% uz = SX.sym('uz');
% controls = [ux;uz]; 
% %nControls = length(controls);
% 
% params=getParams();
% 
% rhs=f_func(states,controls,params);
% f = Function('f', {states,controls}, {rhs}); %nonlinear mapping function f(x,u)
% measurement_rhs = [x;z;theta]; 
% nMeasurements= length(measurement_rhs);
% h = Function('h', {states},{measurement_rhs}); %nonlinear mapping function h(x), @ R^[x,z,theta,x_dot,z_dot,theta_dot] --> R^[x,z,theta]

% mhe=MHEclass(N_MHE,states,controls, Ac,Bc,C,Q,R,M,z0block,x0,dt,CasADiPath,f,h);
% 
% %mhe.xCurrent
% [a,b]=mhe.linearizeDynamics(x0)







if true

    %%%% Run MHE %%%%%%%

    mhe=MHEclass_linear(N_MHE,Ac,Bc,C,Q,R,M,z0block,x0,dt);
    iter_class=0;
    
    %Buffer first horizon of measurements and control inputs
    while ~mhe.isReadyToRun
        newY=Y_noisy(mhe.yBufferCount);
        newU=U_list(mhe.uBufferCount);
        mhe=mhe.bufferInitialData(newY, newU);
    end

    % Run class version of MHE
    xsol2=zeros(nStates,size(Y_noisy,2)-(N_MHE+1));
    for k=1:size(Y_noisy,2)-(N_MHE+1)
        newY=Y_noisy(:,N_MHE+1+k);
        newU=U_list(N_MHE+k,:);
        mhe=mhe.runMHE(newY,newU);
        xsol2(:,k)=mhe.xCurrent;
        
        %run mpc(xk)
        %uk = zopt
        
        iter_class=iter_class+1;

    end
     
     run("plotting.m")
 
end


function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end