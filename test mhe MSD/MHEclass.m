classdef MHEclass
    
    properties
        N_MHE
        nStates
        nControls
        nMeasurements
        z0block
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
        dt
        isReadyToRun
        yBufferCount
        uBufferCount
    end
    
    methods
        function obj = MHEclass(N_MHE, A, B, C,z0block, x0, dt)
            obj.N_MHE = N_MHE;
            obj.A = A;
            obj.B = B;
            obj.z0block=z0block;
            obj.C = C;
            obj.dt = dt;
            obj.x0 = x0;
            obj.nStates= size(A,1);
            obj.nControls = size(B,2);
            obj.nMeasurements=size(C,1);
            
            
            obj.isReadyToRun = false;
            obj.yBufferCount = 0;
            obj.uBufferCount = 0;
            obj = obj.initialize();
        end
        
        function obj=initialize(obj)
            % Setup general stuff, discretizing dynamics and weights etc. 
            obj.A = expm(obj.A * obj.dt);
            obj.B = inv(obj.A) * (expm(obj.A * obj.dt) - eye(size(obj.A))) * obj.B;
            
            obj.Q = diag([0.3, 3]);
            obj.R = diag([0.001]);
            obj.M = diag([0.05, 0.03]);
            
            obj.Q = obj.dt * obj.Q;
            obj.R = obj.dt * obj.R;
            obj.M = obj.dt * obj.M;
            
            obj.P = zeros(obj.nStates+obj.nMeasurements+obj.nControls, 1+(obj.N_MHE+1)+obj.N_MHE); 
            obj.P(1:obj.nStates,1)=obj.x0;

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
        
        
        
        
%         function obj=bufferInitialData(obj, newY, newU)
%             obj.yBufferCount=obj.yBufferCount+1;
%             obj.uBufferCount=obj.uBufferCount+1;
%             
%             
%             %Buffer new measurement
%             if obj.yBufferCount <= obj.N_MHE + 1 && ~isempty(newY)
%                 obj.P(obj.nStates + 1:obj.nStates + obj.nMeasurements, obj.yBufferCount + 1) = newY;
%                 obj.beq(obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.yBufferCount-1)+1 : obj.nStates*obj.N_MHE+obj.nMeasurements*(obj.yBufferCount)) = obj.P(obj.nStates+1:obj.nStates+obj.nMeasurements, obj.yBufferCount+1);
%             end
%             if obj.uBufferCount <= obj.N_MHE && ~isempty(newU)
%                % objP(nStates+nMeasurements+1:nStates+nMeasurements+nControls,1+(N_MHE+1)+k)=U_list(k);
%                 %beq(nStates*k-1 : nStates*(k)) = z0_block + B*P(nStates+nMeasurements+1:nStates+nMeasurements+nControls,(N_MHE+1)+1+k: (N_MHE+1)+1+nControls*k);
%             end
%         end
    end

end