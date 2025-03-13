classdef MHEclass_KF_Update
    
    properties
        N_MHE           %Horizon length #timesteps
        nStates         %Number of states
        nControls       %Number of control inputs
        nMeasurements   %Number of measurements
        Ac              %Continuous time system matrix, A
        Bc              %Continuous time control matrix, B
        A               %Discrete time system matrix, A
        B               %Discrete time control matrix, B
        C               %Measurement matrix, C
        Q               %Weight matrix for process noise, w
        R               %Weight matrix for measurement noise, v
        M               %Weight matrix for arrival cost, x0 - xprior
        P               %Parameters matrix. Stores xprior, measurements and control inputs
        G               %QP cost matrix for quadratic terms
        g               %QP cost vector for linear terms
        Aeq             %QP equality constraint matrix for z, Aeq*z=beq
        beq             %QP equality constraint vector for constant terms
        f               %Nonlinear dynamical model, x_{k+1} = f(x_k, u_k) + w_k
        h               %Nonlinear measurement model y_k = h(x_k, u_k) + v_k
        x0              %Initial condition for states
        xprior          %Prior for states (initial condition or previous state estimate)
        P0              %State covariance
        xCurrent        %Current state estimate, xk
        zOpt            %Solution to QP
        dt              %Sampling time for discretization of continuous time dynamics
        isReadyToRun    %Flag to start running MHE
        yBufferCount    %Counter for buffering of measurements in P
        uBufferCount    %Counter for buffering of control inputs in P
        options         %Quadprog solver options
        
        
    end
    
    methods
        function obj = MHEclass_KF_Update(N_MHE, Ac, Bc, C,Q,R,M, x0,P0, dt,options)
            
            %Assigning arguments to class properties
            obj.N_MHE = N_MHE;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.Q=Q;
            obj.R=R;
            obj.M=M;
            obj.C = C;
            obj.dt = dt;
            obj.x0 = x0;
            obj.P0 = P0;
            obj.nStates= size(Ac,1);
            obj.nControls = size(Bc,2);
            obj.nMeasurements=size(C,1);
            obj.options= options;

           
            % running setup
            obj = obj.setup();
        end        
        
        function obj=setup(obj)
            
            % Discretizing dynamics and initializing P. 
            obj.A = expm(obj.Ac * obj.dt);
            obj.B = inv(obj.Ac) * (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * obj.Bc;
            obj.P = zeros(obj.nStates+obj.nMeasurements+obj.nControls, 1+(obj.N_MHE+1)+obj.N_MHE); 

            % Buffering init
            obj.isReadyToRun = false;
            obj.yBufferCount = 1;
            obj.uBufferCount = 1;
            
            
            %Setup optimization problem, G, g, Aeq, beq etc.
            obj = obj.setupOptimizationProblem(); 

        end
        
        function obj=setupOptimizationProblem(obj)
            
            % Constructing G
            cost_X_block = blkdiag(obj.M,zeros(obj.nStates*obj.N_MHE)); %Arrival cost, for x0
            cost_Q_block = kron(obj.Q,eye(obj.N_MHE)); %Stage costs for process noise, w
            cost_R_block= kron(obj.R,eye(obj.N_MHE+1)); %stage costs for measurement noise, v
            obj.G=blkdiag(cost_X_block,cost_Q_block,cost_R_block);
            
            
            % Constructing g
            g_X=[-2*obj.M*obj.P(1:obj.nStates,1);zeros(obj.nStates*obj.N_MHE,1)]; %Arrival cost. Only linear term is -2*M*xprior*x(0)
            g_W = zeros(obj.nStates * obj.N_MHE, 1);
            g_V = zeros(obj.nMeasurements * (obj.N_MHE + 1), 1);
            obj.g = [g_X; g_W; g_V];
            
            % Constructing Aeq and beq
            rows_Aeq= obj.N_MHE*obj.nStates+(obj.N_MHE+1)*obj.nMeasurements;
            cols_Aeq= (obj.N_MHE+1)*obj.nStates + obj.N_MHE*obj.nStates + (obj.N_MHE+1)*obj.nMeasurements;
            obj.Aeq =zeros(rows_Aeq,cols_Aeq);
            obj.beq = zeros(rows_Aeq, 1);
            
            for k = 0:obj.N_MHE-1
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*k+1 : obj.nStates*(k+1)) = -obj.A;
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*(k+1)+1 : obj.nStates*(k+2)) = eye(obj.nStates);
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*(obj.N_MHE+1) + obj.nStates*k + 1 :obj.nStates*(obj.N_MHE+1) + obj.nStates*(k+1)) = -eye(obj.nStates);
            end
            
            for k=0:obj.N_MHE
                obj.Aeq(obj.nStates*obj.N_MHE+obj.nMeasurements*k+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(k+1), obj.nStates*k+1 : obj.nStates*(k+1))=obj.C;
                obj.Aeq(obj.nStates*obj.N_MHE+obj.nMeasurements*k+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(k+1), obj.nStates*(obj.N_MHE+1)+obj.nStates*(obj.N_MHE)+obj.nMeasurements*k+1 : obj.nStates*(obj.N_MHE+1)+obj.nStates*(obj.N_MHE)+obj.nMeasurements*(k+1)) = eye(obj.nMeasurements);
            end
        end
        
        function obj=bufferInitialData(obj, newY, newU)
            
            %Checking if the horizon is filled
            if obj.yBufferCount > obj.N_MHE + 1 && obj.uBufferCount > obj.N_MHE
                obj=initialGuessPropegation(obj); %Start propagating initial guess over the first horizon
                obj.isReadyToRun = true; % Indicate that the system is ready to start running the MHE
                return
            end
            
            %Buffer new measurement
            if obj.yBufferCount <= obj.N_MHE + 1 && ~isempty(newY)
                obj.P(obj.nStates + 1:obj.nStates + obj.nMeasurements, obj.yBufferCount + 1) = newY;
                obj.beq(obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.yBufferCount-1)+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.yBufferCount)) = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, obj.yBufferCount+1);
                obj.yBufferCount=obj.yBufferCount+1;
            end
           
            %Buffer new control input
            if obj.uBufferCount <= obj.N_MHE && ~isempty(newU)
                obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+obj.uBufferCount)=newU;
                obj.beq(obj.nStates*(obj.uBufferCount-1)+1 : obj.nStates*obj.uBufferCount) =  obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,(obj.N_MHE+1)+1+obj.uBufferCount);
                obj.uBufferCount=obj.uBufferCount+1;
            end
            %^ A biproduct of this solution is that the bufferCounts will
            %always be one more after its finished than the stopping value
            %because after the last iteration it increases once more.
            
        end
        
        function obj=initialGuessPropegation(obj)
            
            x_propagated = obj.x0;

            %The following is not currently in use and is probably useless
                %Integrating x0 over the first horizon before MHE starts
            for i = 1:obj.N_MHE
                %x_propagated = obj.A * x_propagated + obj.B * obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+i);
            end
            obj.xprior=x_propagated;
            obj.P(1:obj.nStates, 1) = x_propagated; %Update P with the propagated initial condition
            g_X(1:obj.nStates)=-2*obj.M*obj.P(1:obj.nStates,1); %Update g in cost function with the propagated initial condition
            obj.g(1:obj.nStates)=g_X(1:obj.nStates);
        end
        
        function obj=runMHE(obj,newY,newU)

            
            
            %Shift measurement window with new measurement
            obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2:obj.N_MHE+2)=[obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 3:obj.N_MHE+2),newY]; 
    
            %Shift control input window with new control input
            obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1:1+(obj.N_MHE+1)+obj.N_MHE)=[obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1+1:1+(obj.N_MHE+1)+obj.N_MHE),newU]; 
            
            %update the C*Xk + Vk = Ymeas,k constraint with new measurement in beq
            obj.beq(obj.nStates*obj.N_MHE+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.N_MHE+1)) = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2:obj.N_MHE+2); 

            %Update dynamics constraint with new control input in beq
            obj.beq(1:obj.nStates*obj.N_MHE)=reshape(obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1:1+(obj.N_MHE+1)+obj.N_MHE),obj.nStates*obj.N_MHE,1);

            %solve opt problem, currently only extracting zOpt
            [obj.zOpt, ~, ~, ~,~] = quadprog(2*obj.G, obj.g,[],[], obj.Aeq, obj.beq, [], [], [], obj.options);
            
            x_zero = obj.zOpt(1:obj.nStates);
            % Extract the control input at the start of the horizon (k-N)
            u_kN = obj.P(obj.nStates+obj.nMeasurements+1:end, 1+(obj.N_MHE+1)+1); % Correct index for u_{k-N}
            
            % Predict x_zero (k-N) to k-N+1 using u_{k-N}
            x_predicted = obj.A * x_zero + obj.B * u_kN; 
            P_predicted = obj.A * obj.P0 * obj.A' + inv(obj.Q); % Q_MHE is inverse covariance
            
            % Extract measurement at k-N+1 (within horizon, not newY)
            y_kN_plus1 = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2); % Index 2 corresponds to k-N+1
            
            % Compute innovation using the correct measurement
            innovation = y_kN_plus1 - obj.C * x_predicted;
            
            % Kalman gain and update
            S_cov_upd = obj.C * P_predicted * obj.C' + inv(obj.R);
            % Regularize S_cov_upd to avoid numerical issues
            S_reg = S_cov_upd + 1e-6 * eye(size(S_cov_upd));
            kalman_gain = P_predicted * obj.C' / S_reg; % Use matrix division for stability
            
            x_corrected = x_predicted + kalman_gain * innovation;
            P_corrected = (eye(obj.nStates) - kalman_gain * obj.C) * P_predicted;
            
            % Regularize P_corrected before inverting
            P_reg = P_corrected + 1e-6 * eye(obj.nStates);
            obj.P0 = P_reg;
            newM = inv(P_reg);
            %obj.M = inv(P_reg); % Update arrival cost weight
            
            % Update xprior for next iteration
            %obj.xprior = x_corrected;
            % x_zero = obj.zOpt(1:obj.nStates);
            % x_predicted = obj.A*x_zero + obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1)
            % P_predicted = obj.A*obj.P0*transpose(obj.A) + inv(obj.Q)
            % innovation = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements,1+obj.N_MHE+1)-obj.C*x_predicted
            % S_cov_upd = obj.C * P_predicted * transpose(obj.C) + inv(obj.R)
            % kalman_gain = P_predicted*transpose(obj.C)*inv(S_cov_upd)
            % x_corrected = x_predicted + kalman_gain * innovation;
            % P_corrected = (eye(size(x_predicted))- kalman_gain * obj.C) * P_predicted;
            % obj.P0=P_corrected;
            % 
            % epsilon = 1e-6;  % or another small value
            % P_reg = P_corrected + epsilon * eye(size(P_corrected));
            % newM = inv(P_reg)           
            % %obj.M=newM;
            % 
            % obj.xprior = x_corrected;
            %obj.xprior = obj.zOpt(obj.nStates + 1 : 2 * obj.nStates);  %Update xprior as the next element after x0 for next iteration
            obj.xprior = x_corrected ;

            %priorError=obj.xprior - x_corrected
            newM=1/2*(newM+transpose(newM));
            obj.M=newM;
            %obj.xprior
            %x_corrected

            obj.xCurrent=obj.zOpt(obj.nStates*obj.N_MHE + 1 : obj.nStates*(obj.N_MHE + 1)); %Extract current state estimate xk as the last element
            % Update arrival cost with new xprior
            obj.P(1:obj.nStates,1) = obj.xprior; 
            g_X(1:obj.nStates) = -2*obj.M*obj.P(1:obj.nStates,1);
            obj.g(1:obj.nStates) = g_X(1:obj.nStates);

            %Update M block in G
            cost_X_block = blkdiag(obj.M,zeros(obj.nStates*obj.N_MHE)); %Arrival cost, for x0
            cost_Q_block = kron(obj.Q,eye(obj.N_MHE)); %Stage costs for process noise, w
            cost_R_block= kron(obj.R,eye(obj.N_MHE+1)); %stage costs for measurement noise, v
            obj.G=blkdiag(cost_X_block,cost_Q_block,cost_R_block);

            
        end
        
        function obj = reset(obj, newx0)
            
            %Reset initial condition, xCurrent and run setup again
            obj.x0 = newx0;           
            obj.xCurrent = [];
            obj = obj.setup();  % Reinitialize the matrices if necessary
        end
        
    end

end