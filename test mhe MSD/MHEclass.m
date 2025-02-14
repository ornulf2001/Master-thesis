classdef MHEclass
    
    properties
        N_MHE
        nStates
        nControls
        nMeasurements
        z0block
        Ac
        Bc
        A
        B
        C
        Q
        R
        M
        P
        G
        g
        Aeq
        beq
        x0
        xprior
        xCurrent
        wCurrent
        vCurrent
        dt
        isReadyToRun
        yBufferCount
        uBufferCount
        lenZ
        nWSR
        
    end
    
    methods
        function obj = MHEclass(N_MHE, Ac, Bc, C,Q,R,M,z0block, x0, dt)
            
            %Assigning arguments to class properties
            obj.N_MHE = N_MHE;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.Q=Q;
            obj.R=R;
            obj.M=M;
            obj.z0block=z0block;
            obj.C = C;
            obj.dt = dt;
            obj.x0 = x0;
            obj.nStates= size(Ac,1);
            obj.nControls = size(Bc,2);
            obj.nMeasurements=size(C,1);
            obj.lenZ=obj.nStates*(obj.N_MHE+1)+obj.nStates*(obj.N_MHE)+obj.nMeasurements*(obj.N_MHE+1);
           
            % Solver options
            obj.nWSR=1000;
            
            % running setup
            obj = obj.setup();
        end
        
        function obj=setup(obj)
            
            % Setup general stuff, discretizing dynamics and weights etc. 
            obj.A = expm(obj.Ac * obj.dt);
            obj.B = inv(obj.Ac) * (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * obj.Bc;
         
            
            obj.P = zeros(obj.nStates+obj.nMeasurements+obj.nControls, 1+(obj.N_MHE+1)+obj.N_MHE); 

            % Buffering init
            obj.isReadyToRun = false;
            obj.yBufferCount = 1;
            obj.uBufferCount = 1;
            
            %Setup opt problem, G, g, Aeq, beq etc.
            obj = obj.setupOptimizationProblem(); 

        end
        
        function obj=setupOptimizationProblem(obj)
            
            % Constructing G
            cost_X_block = blkdiag(obj.M,zeros(obj.nStates*obj.N_MHE));
            cost_Q_block = kron(obj.Q,eye(obj.N_MHE));
            cost_R_block= kron(obj.R,eye(obj.N_MHE+1));
            obj.G=blkdiag(cost_X_block,cost_Q_block,cost_R_block);
            
            
            % Constructing g
            g_X=[-2*obj.M*obj.P(1:obj.nStates,1);zeros(obj.nStates*obj.N_MHE,1)]; %Only linear term is -2*M*x_prior*x(0)
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
                % Indicate that the system is ready to start running the MHE
                obj=initialGuessPropegation(obj);
                obj.isReadyToRun = true;
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
                obj.beq(obj.nStates*obj.uBufferCount-1 : obj.nStates*obj.uBufferCount) = obj.z0block + obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,(obj.N_MHE+1)+1+obj.uBufferCount);
                obj.uBufferCount=obj.uBufferCount+1;
            end
            %^ A biproduct of this solution is that the bufferCount's will
            %always be one more after its finished than the stopping value
            %because after the last iteration it increases once more.
            
        end
        
        function obj=initialGuessPropegation(obj)
            
            x_propagated = obj.x0;
            for i = 1:obj.N_MHE
                x_propagated = obj.A * x_propagated + obj.B * obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+i);
            end
            obj.P(1:obj.nStates, 1) = x_propagated;
            g_X(1:obj.nStates)=-2*obj.M*obj.P(1:obj.nStates,1);
            obj.g(1:obj.nStates)=g_X(1:obj.nStates);
        end
        
        function obj=runMHE(obj,newY,newU)
            
            %solve opt problem, currently only extracting zOpt
            [zOpt, ~, ~, ~,~,~] = qpOASES(2*obj.G, obj.g, obj.Aeq, [],[], obj.beq, obj.beq,[]);
            obj.xprior = zOpt(obj.nStates + 1 : 2 * obj.nStates);  
            obj.xCurrent=zOpt(obj.nStates*obj.N_MHE + 1 : obj.nStates*(obj.N_MHE + 1));
            % Update arrival cost with new xprior
            obj.P(1:obj.nStates,1) = obj.xprior; 
            g_X(1:obj.nStates) = -2*obj.M*obj.P(1:obj.nStates,1);
            obj.g(1:obj.nStates) = g_X(1:obj.nStates);
            
            %Shift measurement window
            obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2:obj.N_MHE+2)=[obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 3:obj.N_MHE+2),newY]; 
    
            %Shift control input window
            obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1:1+(obj.N_MHE+1)+obj.N_MHE)=[obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1+1:1+(obj.N_MHE+1)+obj.N_MHE),newU]; 
            
            %update the C*Xk + Vk = Ymeas,k constraint in beq with new measurement
            obj.beq(obj.nStates*obj.N_MHE+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.N_MHE+1)) = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, 2:obj.N_MHE+2); 

            %Update dynamics constraint with new control input in beq
            obj.beq(1:obj.nStates*obj.N_MHE)=reshape(obj.z0block+obj.B*obj.P(obj.nStates+obj.nMeasurements+1:obj.nStates+obj.nMeasurements+obj.nControls,1+(obj.N_MHE+1)+1:1+(obj.N_MHE+1)+obj.N_MHE),obj.nStates*obj.N_MHE,1);
        end
        
    end

end