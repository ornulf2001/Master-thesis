%Initialization, fill up the first horizon of measurements and control inputs

nWSR=1000; %Max iterations
for k=1:N_MHE+1 %Fill up first horizon of measurements before loop in P and then in beq
    P(nStates+1:nStates+nMeasurements, k+1)=Y_noisy(k);
    beq(nStates*N_MHE+nMeasurements*(k-1)+1 : nStates*N_MHE+nMeasurements*(k)) = P(nStates+1:nStates+nMeasurements, k+1);
end
for k=1:N_MHE %Fill up first horizon-1 of control inputs applied in P and then in beq
    P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+k)=U_list(k);
    beq(nStates*k-1 : nStates*k) = z0_block + B*P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+k);
end
 x_propagated = x0_sim;
for i = 1:N_MHE
    x_propagated = A * x_propagated + B * P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+i);
end
% Update prior in P matrix
P(1:nStates, 1) = x_propagated;
g_X(1:nStates)=-2*M*P(1:nStates,1);
g(1:nStates)=g_X(1:nStates);