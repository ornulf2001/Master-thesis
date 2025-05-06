classdef MHEclass_KF_Update
    
    properties
        N             %Horizon length #timesteps
        nStates           %Number of states
        nControls         %Number of control inputs
        nMeasurements     %Number of measurements
        Ac                %Continuous time system matrix, A
        Bc                %Continuous time control matrix, B
        A                 %Discrete time system matrix, A
        B                 %Discrete time control matrix, B
        C                 %Measurement matrix, C
        Q                 %Weight matrix for process noise, w
        R                 %Weight matrix for measurement noise, v
        M                 %Weight matrix for arrival cost, x0 - xprior
        weightScaling     %Downscaling of Q, R, M matrices used in the filter to avoid ill posed hessian
        P                 %Parameters matrix. Stores xprior, measurements and control inputs
        G                 %QP cost matrix for quadratic terms
        g                 %QP cost vector for linear terms
        Aeq               %QP equality constraint matrix for z, Aeq*z=beq
        beq               %QP equality constraint vector for constant terms
        lb                %Decision variable constraints lb
        ub                % ---||--- ^ ub
        f                 %Nonlinear dynamical model, x_{k+1} = f(x_k, u_k) + w_k
        h                 %Nonlinear measurement model y_k = h(x_k, u_k) + v_k
        x0                %Initial condition for states
        xlp               %Linearization point for nonlinear model    
        xprior            %Prior for states (initial condition or previous state estimate)
        P0                %State covariance
        currentP            %Forward propegated state estimate covariance at k
        xCurrent          %Current state estimate, xk
        vCurrent          %Current measurement noise estimate
        wCurrent          %Current process noise estimate
        currentNIS        %Current NIS
        currentInnovation %Current innovation
        zOpt              %Solution to QP
        zPrev             %Previous solution to QP for warm start
        lastWsetB         %Previous working set - bounds - for warm start
        lastWsetC         %Previous working set - constraints - for warm start
        dt                %Sampling time for discretization of continuous time dynamics
        isReadyToRun      %Flag to start running MHE
        yBufferCount      %Counter for buffering of measurements in P
        uBufferCount      %Counter for buffering of control inputs in P
        options           %Quadprog solver options
        iterations        %Number of iterations run
        
        
    end
    
    methods
        function obj = MHEclass_KF_Update(N, Ac, Bc, C,Q,R,M,weightScaling, x0,xlp,P0, dt,options)
            
            %Assigning arguments to class properties
            obj.N = N;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.Q=Q;
            obj.R=R;
            obj.M=M;
            obj.weightScaling = weightScaling;
            obj.C = C;
            obj.dt = dt;
            obj.xlp = xlp;
            obj.x0 = x0;
            obj.P0 = P0;
            obj.nStates= size(Ac,1);
            obj.nControls = size(Bc,2);
            obj.nMeasurements=size(C,1);
            obj.options= options;
            obj.iterations = 0;

           
            % running setup
            obj = obj.setup();
        end        
        
        function obj=setup(obj)
            
            % Discretizing dynamics and initializing P. 
            %obj.A = expm(obj.Ac * obj.dt);
            %obj.B = (obj.A - eye(size(obj.Ac))) * (obj.Ac \ obj.Bc);
            %obj.B = inv(obj.Ac) * (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * obj.Bc;
            obj.A=eye(size(obj.Ac))+obj.Ac*obj.dt;
            obj.B=obj.Bc*obj.dt;

            obj.P = zeros(obj.nStates+obj.nMeasurements+obj.nControls, 1+(obj.N+1)+obj.N); 
            obj.P(1:obj.nStates,1)=obj.x0;
            obj.P(obj.nStates + 1:obj.nStates + obj.nMeasurements,2:obj.N+2)=(obj.C*obj.x0).*ones(obj.nMeasurements,obj.N+1);
            %obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2) = obj.C*obj.x0;
            obj.zPrev = [obj.x0;zeros(obj.N*obj.nStates,1);zeros(obj.N*obj.nStates,1);zeros((obj.N+1)*obj.nMeasurements,1)];

            nVars = length(obj.zPrev);          % Total number of decision variables
            nConstr = size(obj.Aeq, 1);         % Total number of general constraints

            % On first iteration only:
            obj.lastWsetB = zeros(nVars, 1);    % No bounds active
            obj.lastWsetC = zeros(nConstr, 1);  % No constraints active

            % Buffering init
            obj.isReadyToRun = false;
            obj.yBufferCount = 1;
            obj.uBufferCount = 1;
            
            
            %Setup optimization problem, G, g, Aeq, beq etc.
            obj = obj.setupOptimizationProblem(); 

        end
        
        function obj=setupOptimizationProblem(obj)
            
            % Constructing G
            cost_X_block = blkdiag(obj.weightScaling*obj.M,zeros(obj.nStates*obj.N)); %Arrival cost, for x0
            cost_Q_block = kron(eye(obj.N),obj.weightScaling*obj.Q); %Stage costs for process noise, w
            cost_R_block= kron(eye(obj.N+1),obj.weightScaling*obj.R); %stage costs for measurement noise, v
            obj.G=blkdiag(cost_X_block,cost_Q_block,cost_R_block);
            
            
            % Constructing g
            g_X=[-2*obj.weightScaling*obj.M*obj.P(1:obj.nStates,1);zeros(obj.nStates*obj.N,1)]; %Arrival cost. Only linear term is -2*M*xprior*x(0)
            g_W = zeros(obj.nStates * obj.N, 1);
            g_V = zeros(obj.nMeasurements * (obj.N + 1), 1);
            obj.g = [g_X; g_W; g_V];
            
            % Constructing Aeq and beq
            rows_Aeq= obj.N*obj.nStates+(obj.N+1)*obj.nMeasurements;
            cols_Aeq= (obj.N+1)*obj.nStates + obj.N*obj.nStates + (obj.N+1)*obj.nMeasurements;
            obj.Aeq =zeros(rows_Aeq,cols_Aeq);
            obj.beq = zeros(rows_Aeq, 1);
            
            for k = 0:obj.N-1
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*k+1 : obj.nStates*(k+1)) = -obj.A;
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*(k+1)+1 : obj.nStates*(k+2)) = eye(obj.nStates);
                obj.Aeq(obj.nStates*k+1 : obj.nStates*(k+1), obj.nStates*(obj.N+1) + obj.nStates*k + 1 :obj.nStates*(obj.N+1) + obj.nStates*(k+1)) = -eye(obj.nStates);
            end
            
            for k=0:obj.N
                obj.Aeq(obj.nStates*obj.N+obj.nMeasurements*k+1 : obj.nStates*obj.N+obj.nMeasurements*(k+1), obj.nStates*k+1 : obj.nStates*(k+1))=obj.C;
                obj.Aeq(obj.nStates*obj.N+obj.nMeasurements*k+1 : obj.nStates*obj.N+obj.nMeasurements*(k+1), obj.nStates*(obj.N+1)+obj.nStates*(obj.N)+obj.nMeasurements*k+1 : obj.nStates*(obj.N+1)+obj.nStates*(obj.N)+obj.nMeasurements*(k+1)) = eye(obj.nMeasurements);
            end
            obj.Aeq=sparse(obj.Aeq);
            obj.G=sparse(obj.G);


            nX = obj.nStates;
            zDim = length(obj.zPrev);   % total number of decision variables
            lb = -inf(zDim, 1)         % default lower bound
            ub =  inf(zDim, 1);         % default upper bound
            
            % Constrain z-position (index 3) of all x0...xN to within [z_min, z_max]
            z_min = 0.02;
            z_max = 0.1;
            
            for j = 0:obj.N
                idx = j*nX + 3;  % index for z position in x_k
                lb(idx) = z_min;
                ub(idx) = z_max;
            end
            
            % Optional: limit zdot (index 8)
            zdot_min = -2;
            zdot_max = 2;
            for j = 0:obj.N
                idx = j*nX + 8;
                lb(idx) = zdot_min;
                ub(idx) = zdot_max;
            end
            
            obj.lb = lb;
            obj.ub = ub;
        end
        
        function obj=bufferInitialData(obj, newY, newU) %% This is not in use
            
            %Checking if the horizon is filled
            if obj.yBufferCount > obj.N + 1 && obj.uBufferCount > obj.N
                obj.isReadyToRun = true; % Indicate that the system is ready to start running the MHE
                return
            end
            
            %Buffer new measurement
            if obj.yBufferCount <= obj.N + 1 && ~isempty(newY)
                obj.P(obj.nStates + 1:obj.nStates + obj.nMeasurements, obj.yBufferCount + 1) = newY;
                obj.beq(obj.nStates*obj.N+obj.nMeasurements*(obj.yBufferCount-1)+1 : obj.nStates*obj.N+obj.nMeasurements*(obj.yBufferCount)) = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, obj.yBufferCount+1);
                obj.yBufferCount=obj.yBufferCount+1;
            end
           
            %Buffer new control input
            if obj.uBufferCount <= obj.N && ~isempty(newU)
                obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N+1)+obj.uBufferCount)=newU;
                obj.beq(obj.nStates*(obj.uBufferCount-1)+1 : obj.nStates*obj.uBufferCount) =  obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,(obj.N+1)+1+obj.uBufferCount);
                obj.uBufferCount=obj.uBufferCount+1;
            end
            %^ A biproduct of this solution is that the bufferCounts will
            %always be one more after its finished than the stopping value
            %because after the last iteration it increases once more.
            
        end
        
        function obj = initialGuessPropegation(obj) %% This is not in use

            % Start from your known initial state x0
            x_propagated = obj.x0;
        
            % Example: fill the 'measurement guess' for the first sample
            %   Suppose columns 2..(N+2) store the horizon's measured outputs
            obj.P(obj.nStates+1 : obj.nStates+obj.nMeasurements, 2:1+obj.N+1) = obj.C * x_propagated;
        
            for i = 1 : obj.N
                % 1) Choose a nominal guess for the control (often zero):
                c_i = zeros(obj.nControls,1);
        
                % 2) Store this nominal control in the correct column
                %    e.g. columns (N+2)..(2*N+1) might hold controls:
                colControl = (obj.N+1) + 1 + (i);
                obj.P(obj.nStates+obj.nMeasurements+1 : end, colControl) = c_i;
        
                % 3) Forward-simulate to get x_{i+1}
                x_propagated = obj.A * x_propagated + obj.B * c_i;
                if i==1
                    obj.P(1:obj.nStates,1)=x_propagated;
                end
        
                % 4) Store a measurement guess for time i+1
                colMeas = 2 + i;  % next measurement column
                obj.P(obj.nStates+1 : obj.nStates+obj.nMeasurements, colMeas) = obj.C * x_propagated;
            end
        end
        
        function obj=runMHE(obj,newY,newU)
            obj.iterations = obj.iterations+1;
            

            % Row and column indices for control-block and measurement-block in P
            ctrl_start = obj.nStates + obj.nMeasurements + 1;
            ctrl_end   = ctrl_start + obj.nControls - 1;
            col_ctrl_start = 1 + (obj.N+1) + 1;
            col_ctrl_end   = 1 + (obj.N+1) + obj.N;
            meas_start = obj.nStates + 1;
            meas_end   = obj.nStates + obj.nMeasurements;
            col_meas_start = 2;
            col_meas_end = 1 + (obj.N+1);  % which equals N + 2

            if obj.N == 1
                % For a one-step horizon, just assign newU to the only control column.
                obj.P(ctrl_start:ctrl_end, col_ctrl_start) = newU;
            else
                % For N > 1, shift left by one column and insert the new control.
                obj.P(ctrl_start:ctrl_end, col_ctrl_start:col_ctrl_end) = ...
                    [ obj.P(ctrl_start:ctrl_end, col_ctrl_start+1:col_ctrl_end), newU ];
            end

            obj.P(meas_start:meas_end, col_meas_start:col_meas_end) = ...
                [ obj.P(meas_start:meas_end, col_meas_start+1:col_meas_end), newY ];

            %disp("P after shift")
            %disp(obj.P)
            %Shift measurement window with new measurement
            %obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2:obj.N+2)=[obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 3:obj.N+2),newY]; 
    
            %Shift control input window with new control input
            %obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N+1)+1:1+(obj.N+1)+obj.N)=[obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N+1)+1+1:1+(obj.N+1)+obj.N),newU]; 
            
            %update the C*Xk + Vk = Ymeas,k constraint with new measurement in beq
            obj.beq(obj.nStates*obj.N+1 : obj.nStates*obj.N+obj.nMeasurements*(obj.N+1)) = reshape( obj.P(meas_start:meas_end, col_meas_start:col_meas_end), [], 1 );

            %Update dynamics constraint with new control input in beq
            %obj.beq(1:obj.nStates*obj.N)=reshape(obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N+1)+1:1+(obj.N+1)+obj.N),obj.nStates*obj.N,1);
            obj.beq(1:obj.nStates*obj.N) = reshape(...
                obj.B * obj.P(ctrl_start:ctrl_end, col_ctrl_start:col_ctrl_end), ...
                    obj.nStates*obj.N, 1);


            obj = obj.computeNIS(); %Find current NIS
           
            %solve opt problem, currently only extracting zOpt
            [obj.zOpt, ~, exitflag, ~,~] = quadprog(2*obj.G, obj.g,[],[], obj.Aeq, obj.beq, [], [], obj.zPrev, obj.options);  
            %disp(exitflag)
            %auxInput=qpOASES_auxInput('x0',obj.zPrev,...
            %                           'guessedWorkingSetB', obj.lastWsetB, ...
            %                           'guessedWorkingSetC', obj.lastWsetC);
            %[obj.zOpt, ~, exitflag, iter, ~,auxOutput] = qpOASES(2*obj.G, obj.g, obj.Aeq, [], [], obj.beq, obj.beq,obj.options,auxInput);
            obj.zPrev = obj.zOpt;
            %obj.lastWsetB = auxOutput.workingSetB;
            %obj.lastWsetC = auxOutput.workingSetC;

            %disp("exitflag: "+ string(exitflag)+", iter: "+ string(iter))
            

            x_zero = obj.zOpt(1:obj.nStates); %x0 estimate | MHE gives x^= [x0^,x1^,...xk^]

            % Extract the control input at the start of the horizon (k-N)
            u_zero = obj.P(obj.nStates+obj.nMeasurements+1:end, 1+(obj.N+1)+1); % Correct index for u_{k-N}
            
            % Predict x_zero (k-N) to k-N+1 using u_{k-N}
            x_predicted = obj.A * x_zero + obj.B * u_zero; 
            P_predicted = obj.A * obj.P0 * transpose(obj.A) + inv(obj.Q); % Q_MHE is inverse covariance
            
            % Extract measurement at k-N+1 
            y_k_Nplus1 = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 3); % Index 3 corresponds to k-N+1
            
            % Compute innovation using the correct measurement
            innovation = y_k_Nplus1 - obj.C * x_predicted;
            
            % Kalman gain and update
            S_cov_upd = obj.C * P_predicted * transpose(obj.C) + inv(obj.R);
            % Regularize S_cov_upd to avoid numerical issues
            S_reg = S_cov_upd+ 1e-6 * eye(size(S_cov_upd));

            kalman_gain = P_predicted * transpose(obj.C)/S_reg; % Use matrix division for stability
            
            x_corrected = x_predicted + kalman_gain * innovation;
            P_corrected = (eye(obj.nStates) - kalman_gain * obj.C) * P_predicted * transpose((eye(obj.nStates)-kalman_gain*obj.C)) + kalman_gain*inv(obj.R)*transpose(kalman_gain);
            
            % Regularize P_corrected before inverting
            P_reg = P_corrected+ 1e-6 * eye(obj.nStates);
            obj.P0 = P_reg; %Update arrival cost covariance initial guess
            newM = obj.weightScaling*inv(P_reg); %Calculate new arrival cost weight
          
            obj.xprior = obj.zOpt(obj.nStates + 1 : 2 * obj.nStates);  %Update xprior in arrival cost as x1^ from MHE
            %obj.xprior = x_corrected; %Could instead update xprior as the x_corrected (propagated KF style from x0^)
            obj.M=1/2*(newM+transpose(newM));

            %Only update M ever 2 iteration, this logic could probably be improved
            % if mod(obj.iterations,1)==0
            %     newM=1/2*(newM+transpose(newM));
            %     obj.M=newM;
            % end
            

            obj.xCurrent=obj.zOpt(obj.nStates*obj.N + 1 : obj.nStates*(obj.N + 1)); %Extract current state estimate xk^
            obj.vCurrent=obj.zOpt(obj.nStates*(obj.N+1) + obj.nStates*obj.N+obj.nMeasurements*obj.N+1: obj.nStates*(obj.N+1) + obj.nStates*obj.N+obj.nMeasurements*(obj.N+1));
            obj.wCurrent=obj.zOpt(obj.nStates*(obj.N+1)+obj.nStates*obj.N+1:obj.nStates*(obj.N+1)+obj.nStates*(obj.N+1) );
            % Update arrival cost with new xprior
            obj.P(1:obj.nStates,1) = obj.xprior; 

            %Update G and g in cost function
            obj.g(1:obj.nStates) = -2*obj.M*obj.P(1:obj.nStates,1);
            obj.G(1:obj.nStates, 1:obj.nStates) = obj.M;

            

            
        end
        
        function obj = reset(obj, newx0) %% This is not in use
            
            %Reset initial condition, xCurrent and run setup again
            obj.x0 = newx0;           
            obj.xCurrent = [];
            obj = obj.setup();  % Reinitialize the matrices if necessary
        end

        function obj = computeNIS(obj) 

            % Initialize propagation with the arrival state.
            x_prop = obj.P(1:obj.nStates,1); 
            P_prop = obj.P0;                  
            
            for i = 1:obj.N
                
                ctrl_col = 1 + (obj.N+1) + i;
                u_i = obj.P(obj.nStates+obj.nMeasurements+1 : obj.nStates+obj.nMeasurements+obj.nControls, ctrl_col);
                
                % Prediction step
                x_pred = obj.A * x_prop + obj.B * u_i;
                P_pred = obj.A * P_prop * obj.A' + inv(obj.Q);
                
                % Update step
                meas_col = 1+i+1 ;
                y_i = obj.P(obj.nStates+1 : obj.nStates+obj.nMeasurements, meas_col);
                
                innovation = y_i - obj.C * x_pred;
                S = obj.C * P_pred * obj.C' + inv(obj.R);
                W = P_pred * obj.C' / S;
                
                x_upd = x_pred + W * innovation;
                P_upd = (eye(obj.nStates) - W * obj.C) * P_pred * (eye(obj.nStates) - W * obj.C)' + W * inv(obj.R) * W';
                
                x_prop = x_upd;
                P_prop = P_upd;
            end
            
            % After this loop, x_prop and P_prop represent the predicted state and covariance at time k.
            obj.currentInnovation = innovation;
            obj.currentNIS = innovation' / S * innovation;
            obj.currentP = P_prop;
            
        end
    end
end