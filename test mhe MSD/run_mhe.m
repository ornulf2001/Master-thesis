for k=1:size(Y_noisy,2)-(N_MHE+1)
[zOpt, fval, exitflag, iterations,lambda,auxOutput] = qpOASES(2*G, g, Aeq, -100000000*ones(length(z),1), 100000000*ones(length(z),1), beq, beq,nWSR); %Solve
xCurrent = zOpt(nStates*N_MHE + 1 : nStates*(N_MHE + 1));  % Extract and store X_N
xsol(:,k) = xCurrent;
zs=[zs,zOpt];
P(1:nStates, 1) = xCurrent; % Update xprior for next iteration

%Shift measurement window
P(nStates+1:nStates+nMeasurements, 2:N_MHE+2)=[P(nStates+1:nStates+nMeasurements, 3:N_MHE+2),Y_noisy(N_MHE+1+k)]; 

%Shift control input window
P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+1:1+(N_MHE+1)+N_MHE)=[P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+1+1:1+(N_MHE+1)+N_MHE),U_list(k+N_MHE)]; 

%update the C*Xk + Vk = Ymeas,k constraint in beq with new measurement window
beq(nStates*N_MHE+1 : nStates*N_MHE+nMeasurements*(N_MHE+1)) = P(nStates+1:nStates+nMeasurements, 2:N_MHE+2); 

%Update dynamics constraint with new control input window in beq
beq(1:nStates*N_MHE)=reshape(z0_block+B*P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+1:1+(N_MHE+1)+N_MHE),nStates*N_MHE,1);

% Update arrival cost with new xprior
g_X(1:nStates)=-2*M*P(1:nStates,1);
g(1:nStates)=g_X(1:nStates);
iter_reg=iter_reg+1;

end