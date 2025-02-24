function [Ac,Bc,C,D,T_sim,U_list,X_sim,Y_sim,Y_noisy]=msd_sim(dt,x0_sim,T_END)
    [Ac,Bc,C,D] = f_func_msd(1,1,0.0001); 
    sys = ss(Ac,Bc,C,D);
    T_sim = 0:dt:T_END-dt; 
    
    U_list = 0.5*sin(T_sim)';%[0:dt:10]'; % Constant force input
    X_sim = zeros(length(x0_sim), length(T_sim)); % State storage
    X_sim(:,1) = x0_sim; % Initial state
    A=expm(Ac*dt);
    B=inv(Ac) * (A - eye(size(Ac))) * Bc;
    % Euler integration for state update
    for k = 1:length(T_sim)-1
        X_sim(:, k+1)=A*X_sim(:, k)+B*U_list(k);
    end
    
    %Synthesizing noisy measurements
    Y_sim=[1,0]*X_sim;
    Y_noisy=Y_sim+0.05*(randn(size(Y_sim)));
end
