clc,clear

N=1:1:3;
dt=0.0001:0.0005:0.001;

[NGrid, dtGrid]=meshgrid(N,dt);
numSimulations = numel(NGrid);
disp(length(N))
    disp(length(dt))
for i=1:numSimulations
    currentN=NGrid(i);
    currentDt=dtGrid(i);
    
    
    disp([currentN,currentDt])
end
    