function [Ac,Bc,z0block,C,D,T_sim,U_list,X_sim,Y_sim,Y_noisy]=sim_2d(dt,x0_sim,T_END,params)
    
    z0=0.0365;
    [Ac,Bc,z0block]=f_func_lin(z0,params);
    C = [eye(3),zeros(3,3)];
    D=0;
    %z0block=zeros(size(Ac,1),1);

    T_sim = 0:dt:T_END-dt; 
    U_list = zeros(size(Bc,2),length(T_sim));%[0:dt:10]'; % Constant force input
    X_sim = zeros(length(x0_sim), length(T_sim)); % State storage
    X_sim(:,1) = x0_sim; % Initial state
    A=expm(Ac*dt);
    B=inv(Ac) * (A - eye(size(Ac))) * Bc;

    xref=[0;0];
    Qk = diag([3,3,6,1,10,10]);
    Rk = diag([10,10]);
    [K,~,~]  = dlqr(A,B,Qk,Rk);
    %F = inv(C * ((Bc*K - Ac) \ Bc));



    % Euler integration for state update
    for k = 1:length(T_sim)-1
        uk=-K*X_sim(:,k); %+ F*xref;
        X_sim(:, k+1)=A*X_sim(:, k)+B*uk;%+z0block;
        U_list(:,k)=uk;
    end
    
    %Synthesizing noisy measurements
    Y_sim=C*X_sim;
    Y_noisy=Y_sim+0.0001*(randn(size(Y_sim)));
    %plot(Y_sim(2,:))
    %ylim([-0.012,0.012])
end

