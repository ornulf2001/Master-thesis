clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*



%%%% Setup %%%%%%%%%


%Dynamical model
z0=0.0365;
params=getParams();
[Ac,Bc,z0block]=f_func_lin(z0,params);

%Measurement model
C=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0];

% Sizes
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);

%
N_MHE=5;
dt=0.005;

R=diag([0.001, 0.001, 0.001]);
Q=diag([1,1,1,1,1,1]);
M=diag([1,1,1,1,1,1]);
x0=[0;0.04;0;0;0;0];






if true

    %%%% Run MHE %%%%%%%

    mhe=MHEclass(N_MHE,Ac,Bc,C,Q,R,M,z0block,x0,dt);
    %mhe=mhe.startLogging(true);
    iter_class=0;
    
    test_Aeq=mhe.Aeq;
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
