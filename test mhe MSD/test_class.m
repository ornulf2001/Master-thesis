clc,clear
run("msd_sim.m")
N_MHE = 10;
z0_block=zeros(size(Ac,1),1);

mhe=MHEclass(N_MHE,Ac,Bc,C,z0_block,x0_sim,dt);

classAeq=mhe.Aeq;

