
addpath ("C:\Users\ornul\AppData\Roaming\MathWorks\CasADi")
import casadi.*

run("mpc_setup.m")

tic
while(norm((x0-xs),2)> 1e-2 && mpcIter < maxIter)
    args.p = [x0;xs]; % Set the values of the parameters vector to the init and setpoint values for the opt variables
    args.x0 = [reshape(X0', 3*(N+1), 1); reshape(u0', 2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
    u = reshape (full(sol.x(3*(N+1)+1:end))',2,N)'; %Get full u from solution
    u_cl(mpcIter+1,:) = u(1,:); %extract and use only first element in control sequence
    
    xx1(:,1:3,mpcIter+1) = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; %Get solution trajectory
    t(mpcIter+1) = t0;
    
    %Apply the control and shift the solution
    [t0, x0, u0, con0] = shift(T, t0, x0, u, f, con_cov);
    xx(:,mpcIter+2) = x0;
    con(mpcIter+1,:)=con0';
    
    X0 = reshape(full(sol.x(1:3*(N+1)))', 3, N+1)'; %Get solution trajectory
    
    %Shift trajectory to initialize the next step
    X0 = [X0(2:end,:); X0(end,:)];
    mpcIter = mpcIter + 1;
end
mpc_iter_time=toc/mpcIter
t_stop_all=find(t==max(t),1);
t_stop=t_stop_all(end);
warning('off');
save("mpc_workspace")

