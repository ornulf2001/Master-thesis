


classdef MPCclass

    properties

        N_MPC
        nStates
        nControls
        zBlock
        Ac
        Bc
        A
        B
        Q
        QN
        R
        G
        Aeq
        beq
        f
        X0
        dt
        lb
        ub
        options
    end

    methods
        function obj = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q, R, QN, nStates, nControls,options)
            obj.N_MPC = N_MPC;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.X0 = X0;
            obj.dt = dt;
            obj.nStates = nStates;
            obj.nControls = nControls;
            obj.Q = Q;
            obj.R = R;
            obj.QN = QN;
            obj.lb = lb;
            obj.ub = ub;
            obj.options = options;
            obj = obj.initialize();
            
        end

        % Function for discretezise A and B matrices
        function obj = initialize(obj)
            obj.A = expm(obj.Ac * obj.dt);
            %obj.B = (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * (obj.Ac \ obj.Bc);
            obj.B = inv(obj.Ac) * (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * obj.Bc;

            % Compute discrete-time B matrix using integral of exp(Ac*t)*Bc from 0 to dt
            %M = expm([-obj.Ac, obj.Bc; zeros(obj.nControls, obj.nStates + obj.nControls)] * obj.dt);
            %Bd = M(1:obj.nStates, obj.nStates+1:end);
            %obj.B = Bd;

            obj = obj.setupOptimizationProblem();

        end

        
        function obj = setupOptimizationProblem(obj)

            % Setting up G and f matrices
            costQ = blkdiag(kron(eye(obj.N_MPC), obj.Q), obj.QN);
            costR = kron(eye(obj.N_MPC), obj.R);
            obj.G = blkdiag(costQ, costR);
            obj.f = zeros(size(obj.G, 1), 1);

            
            % Constructing the Aeq matricx
            nStateBlocks   = obj.N_MPC+1;  % Only x1...x_N_MPC are decision variables.
            nControlBlocks = obj.N_MPC;  % Each step has one control move.
            
            rowAeq = nStateBlocks * obj.nStates;
            colAeq = nStateBlocks * obj.nStates + nControlBlocks * obj.nControls;
            obj.Aeq = zeros(rowAeq, colAeq);
            obj.beq = zeros(rowAeq, 1);
            obj.beq(1:obj.nStates) = obj.X0;
            
             for k = 1:obj.N_MPC+1
                % Row indices for the k-th constraint (each constraint has obj.nStates rows)
                rowIdx = (k-1)*obj.nStates + 1 : k*obj.nStates;
                
                
                colIdx_state = (k-1)*obj.nStates + 1 : k*obj.nStates;
                obj.Aeq(rowIdx, colIdx_state) = eye(obj.nStates);
                
                
                % This will appear “below the diagonal” in the state block.
                if k > 1
                    colIdx_prevState = (k-2)*obj.nStates + 1 : (k-1)*obj.nStates;
                    obj.Aeq(rowIdx, colIdx_prevState) = -obj.A;
                end
                
               
                 if k > 1
                  
                    % Ensure -B is placed only in the correct positions
                    colIdx_control = nStateBlocks*obj.nStates + (k-2)*obj.nControls + 1 : nStateBlocks*obj.nStates + (k-1)*obj.nControls;
                    obj.Aeq(rowIdx, colIdx_control) = -obj.B;
                end
             end
        end

        function [Xopt, Uopt] = runMPC(obj,xCurrent)
        %obj.beq(1:obj.nStates) = xCurrent;
        %b=obj.beq(1:obj.nStates)
         obj.beq(1:obj.nStates) = xCurrent;
        [Zopt, fval, exitflag, output, lambda] = quadprog(2*obj.G, obj.f, [],[],obj.Aeq, obj.beq,obj.lb,obj.ub, [], obj.options);
        
        U0 = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls);
        Uopt = U0(1:obj.nControls);

        X1 = Zopt(1:(obj.N_MPC+1)*obj.nStates);
        xNext=X1(obj.nStates+1:2*obj.nStates);


        Xopt = xNext;
        end



        function plotResults(obj, Xopt, Uopt)

            % Plotting the results
            t = 0:obj.NT;
            length(t)
            
            figure(1);
            subplot(3, 1, 1);
            stairs(t, Xopt(1, :), 'k', 'LineWidth', 0.5);
            grid on;
            legend('$x$', 'Interpreter', 'latex')
            title('States');
            ylabel('$x$', 'Interpreter', 'latex')
            subplot(3, 1, 2);
            stairs(t, Xopt(2, :), 'b',  'LineWidth', 0.5);
            grid on;
            legend('$z$', 'Interpreter', 'latex')
            xlabel('$Time [s]$', 'Interpreter', 'latex')
            ylabel('$z$', 'Interpreter', 'latex')
            subplot(3, 1, 3);
            stairs(t, Xopt(3, :), 'g',  'LineWidth', 0.5);
            grid on;
            legend('$theta$', 'Interpreter', 'latex')
            xlabel('$Time [s]$', 'Interpreter', 'latex')
            ylabel('$theta$', 'Interpreter', 'latex')


            figure(2);
            subplot(2, 1, 1)
            stairs(t(1:end-1), Uopt(1,:), 'r.-',  'LineWidth', 0.5);
            grid on;
            legend('$ux$', 'Interpreter', 'latex');
            ylabel('$ux$', 'Interpreter', 'latex');
            subplot(2, 1, 2)
            stairs(t(1:end-1), Uopt(2,:), 'r.-',  'LineWidth', 0.5);
            grid on;
            legend('$uz$', 'Interpreter', 'latex');
            ylabel('$uz$', 'Interpreter', 'latex');


            figure(3)
            plot(Xopt(1, :), Xopt(2, :));
            grid on;
            ylabel('$x$', 'Interpreter', 'latex');
            xlabel('$z$',  'Interpreter', 'latex');
            title('PhasePortrait of x and z')
        end
    end


end

