
clc,clear
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")                   % Add path to CasADi
addpath ("C:\Users\ornul\Desktop\Kyb master\MASTER\qpOASES\interfaces\matlab")% Add path to qpOAses
addpath (genpath("simulator"))
import casadi.*


%Run msd simulation
x0_sim=[0.1;0];

N=1:1:10;
dt= 0.005;%0.0005:0.001:0.01;
T_END=10;
[NGrid, dtGrid]=meshgrid(N,dt);
nSim = numel(NGrid);

nStates=2;
sims=zeros(nStates,length(0:min(dt):T_END),nSim);
results=zeros(nStates,length(0:min(dt):T_END)-max(N),nSim);
RMSEx=zeros(nSim,1);

for i=1:nSim
    tic
    currentN=NGrid(i);
    currentDt=dtGrid(i);
    disp("Running simulation with N = "+num2str(currentN)+" and dt = "+num2str(currentDt))
    [Ac,Bc,C,D,T_sim,U_list,X_sim,Y_sim,Y_noisy]=msd_sim(currentDt,x0_sim,T_END);
    x0=x0_sim+[0.2;-0.2];
    xsol2=MHE(Ac,Bc,C,currentN,x0,currentDt,Y_noisy,U_list);

    results(1:nStates,1:size(xsol2,2),i)=xsol2;
    sims(1:nStates,1:length(T_sim),i)=X_sim;

    RMSEx(i)=sqrt(mean(xsol2(1,:)-X_sim(1,(currentN+1)+1:size(X_sim,2)).^2));
    toc
    figure(i)
    plot(xsol2(1,:)); hold on
    plot(X_sim(1,currentN+1:size(X_sim,2)-1))
    title("N = "+currentN)
end
disp("Finished all simulations!")
plot(RMSEx)
%RMSEx = reshape(RMSEx, size(NGrid));

% figure(nSim+1)
% surf(NGrid,dtGrid,RMSEx)
% xlabel('Horizon Length N');
% ylabel('Time Step dt');
% zlabel('RMSE');
% title('Surface Plot of RMSE vs. N and dt');
% colorbar;

function xsol=MHE(Ac,Bc,C,N_MHE,x0,dt,Y_noisy,U_list)

    nStates=size(Ac,1);
    z0block=zeros(nStates,1);
    
    
    % Weights and horizon length
    R=0.0005;
    Q=diag([1,1]);
    M=diag([1,1]);


    %%%%%%%%%%% Run MHE %%%%%%%%%%%%%%

    mhe=MHEclass(N_MHE,Ac,Bc,C,Q,R,M,z0block,x0,dt);
    iter_class=0;
    xsol=[];
    
    %Buffer first horizon of measurements and control inputs
    while ~mhe.isReadyToRun
        newY=Y_noisy(mhe.yBufferCount);
        newU=U_list(mhe.uBufferCount);
        mhe=mhe.bufferInitialData(newY, newU);
    end

    % Run MHE
    for k=1:size(Y_noisy,2)-(N_MHE+1)
        newY=Y_noisy(:,N_MHE+1+k);
        newU=U_list(N_MHE+k,:);
        mhe=mhe.runMHE(newY,newU);
        xsol=[xsol,mhe.xCurrent];
        
        iter_class=iter_class+1;
    end
    mhe=mhe.reset(x0);
end
