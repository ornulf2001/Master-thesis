clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*


%Run msd simulation
dt=0.003;
T_END=1;
x0_sim=[0.01;0.04;0.03;0;0;0];
params=getParams();
[Ac,Bc,z0block,C,D,T_sim,U_list,X_sim,Y_sim,Y_noisy]=sim_2d(dt,x0_sim,T_END,params);
z0block=zeros(size(Ac,1),1);
%%%%%%%% Setup %%%%%%%%%

% Sizes
nStates=size(Ac,1);
nControls = size(Bc,2);
nMeasurements = size(C,1);


% Weights and 
N_MHE=10;
R=0.005*diag([1,2,1]);
Q=diag([1,3,1,1,1,1]);
M=diag([3,3,3,1,1,1]);

if true

    %%%%%%%%% Run MHE %%%%%%%%%%

    mhe=MHEclass(N_MHE,Ac,Bc,C,Q,R,M,z0block,x0_sim,dt);
    iter_class=0;
    
    %Buffer first horizon of measurements and control inputs
    while ~mhe.isReadyToRun
        newY=Y_noisy(:,mhe.yBufferCount);
        newU=U_list(:,mhe.uBufferCount);
        mhe=mhe.bufferInitialData(newY, newU);
    end
    Aeq_test=mhe.Aeq;
    beq_test=mhe.beq;
    G_test=mhe.G;
    g_test=mhe.g;
    P_test=mhe.P;
    a=[];
    % Run class version of MHE
    xsol2=zeros(nStates,size(Y_noisy,2)-(N_MHE+1));
    for k=1:size(Y_noisy,2)-(N_MHE+1)
        newY=Y_noisy(:,N_MHE+1+k);
        newU=U_list(:,N_MHE+k);
        mhe=mhe.runMHE(newY,newU);
        xsol2(:,k)=mhe.xCurrent;
        %a=[a,mhe.P];
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
