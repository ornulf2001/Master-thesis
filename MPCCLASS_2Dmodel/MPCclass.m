


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
        xRef
        uRef
        NT
        options
    end

    methods
        function obj = MPCclass(N_MPC, Ac, Bc, xRef, zBlock, X0, dt, lb, ub, Q, R, QN, NT, nStates, nControls)
            obj.N_MPC = N_MPC;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.zBlock = zBlock;
            % obj.xRef = xRef;
            obj.X0 = X0;
            obj.dt = dt;
            obj.nStates = nStates;
            obj.nControls = nControls;
            obj.Q = Q;
            obj.R = R;
            obj.QN = QN;
            obj.lb = lb;
            obj.ub = ub;
            obj.NT = NT;
            obj.xRef = xRef;
            obj = obj.initialize();
        end

        function obj = initialize(obj)
            obj.A = expm(obj.Ac * obj.dt);
            % obj.B = inv(obj.Ac) * (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * obj.Bc;
            obj.B = (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * (obj.Ac \ obj.Bc);


            % obj.A = expm(obj.Ac * obj.dt);
            % M = expm([obj.Ac, obj.Bc; zeros(obj.nControls, obj.nStates + obj.nControls)] * obj.dt);
            % obj.B = M(1:obj.nStates, obj.nStates+1:end);

            % obj.uRef = obj.B \ ((eye(obj.nStates) - obj.A) * obj.xRef);
            obj.uRef = obj.B \ ((eye(obj.nStates) - obj.A)*obj.xRef);

            obj = obj.setupOptimizationProblem();

        end

        function obj = setupOptimizationProblem(obj)
            costQ = blkdiag(kron(eye(obj.N_MPC), obj.Q), obj.QN);
            costR = kron(eye(obj.N_MPC+1), obj.R);
            obj.G = blkdiag(costQ, costR);
            obj.G = obj.G + 1e-6 * eye(size(obj.G));
            

            obj.f = zeros(size(obj.G, 1), 1);

            
            % 
            % stateTerms = repmat(-2 * obj.Q *obj.xRef, obj.N_MPC, 1);
            % termTerm = -2 * obj.QN * obj.xRef;
            % controlTerms = repmat(-2 * obj.R * obj.uRef, obj.N_MPC+1, 1);
            % obj.f = [stateTerms; termTerm; controlTerms];

           % Assumptions:
            %  - obj.nStates: number of states
            %  - obj.nControls: number of controls
            %  - obj.N_MPC: prediction horizon (number of dynamic constraints)
            %
            % Decision variable ordering: 
            %    [x1; x2; ...; x_{N_MPC}; u0; u1; ...; u_{N_MPC-1}]
            %  (x0 is assumed given)
            
            nStateBlocks   = obj.N_MPC+1;  % Only x1...x_N_MPC are decision variables.
            nControlBlocks = obj.N_MPC+1;  % Each step has one control move.
            
            rowAeq = nStateBlocks * obj.nStates;
            colAeq = nStateBlocks * obj.nStates + nControlBlocks * obj.nControls;
            obj.Aeq = zeros(rowAeq, colAeq);
            obj.beq = zeros(rowAeq, 1);
            obj.beq(1:obj.nStates) = obj.A*obj.X0;
            
            for k = 1:obj.N_MPC+1
                % Row indices for the k-th constraint (each constraint has obj.nStates rows)
                rowIdx = (k-1)*obj.nStates + 1 : k*obj.nStates;
                
                % --- Left Block: States ---
                % Place the identity in the k-th state block (this is our "diagonal")
                % Since the state decision variables are ordered as x1, x2, ...,
                % the k-th state variable x_k occupies columns:
                colIdx_state = (k-1)*obj.nStates + 1 : k*obj.nStates;
                obj.Aeq(rowIdx, colIdx_state) = eye(obj.nStates);
                
                % For k > 1, place -A in the block corresponding to the previous state x_{k-1}
                % This will appear “below the diagonal” in the state block.
                if k > 1
                    colIdx_prevState = (k-2)*obj.nStates + 1 : (k-1)*obj.nStates;
                    obj.Aeq(rowIdx, colIdx_prevState) = -obj.A;
                end
                
                % --- Right Block: Controls ---
                % The control blocks follow all the state blocks.
                % For the k-th constraint, the control used is u_{k-1}. Its block is:
                % colIdx_control = nStateBlocks*obj.nStates + (k-1)*obj.nControls + 1 : nStateBlocks*obj.nStates + k*obj.nControls;
                % obj.Aeq(rowIdx, colIdx_control) = -obj.B;
                 if k > 1
                    colIdx_prevState = (k-2)*obj.nStates + 1 : (k-1)*obj.nStates;
                    obj.Aeq(rowIdx, colIdx_prevState) = -obj.A;
            
                    % --- Right Block: Controls (Lower Diagonal) ---
                    % Ensure -B is placed only in the correct positions
                    colIdx_control = nStateBlocks*obj.nStates + (k-2)*obj.nControls + 1 : nStateBlocks*obj.nStates + (k-1)*obj.nControls;
                    obj.Aeq(rowIdx, colIdx_control) = -obj.B;
                end
            end



            % obj.beq = zeros(totalStates, 1);
            % obj.beq(1:obj.nStates) = obj.X0 - obj.xRef;

        end

        function [Xhis, Uopt] = runMPC(obj)
            obj.options = optimset('Display','on', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');

            Xhis = zeros(obj.nStates, obj.N_MPC+1);
            Uopt = zeros(obj.nControls, obj.N_MPC);
            Xhis(:, 1) = obj.X0;
            % 
            % totalStates = (obj.N_MPC + 1) * obj.nStates;
            % totalControls = obj.N_MPC * obj.nControls;

            for k = 1:obj.NT


                % [Zopt, fval, exitflag, output, lambda] = quadprog(obj.G, [], [], [], obj.Aeq, obj.beq, obj.lb, obj.ub, [], obj.options);
                % disp(Zopt)

                % 
                obj.beq(1:obj.nStates) = Xhis(:, k);

                [Zopt, ~, ~, ~, ~, ~] = qpOASES(2*obj.G, obj.f, obj.Aeq, obj.lb, obj.ub, obj.beq, obj.beq, []);

                % U = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls);
                % % disp(U)
                % % disp(dU)
                % % U = dU + obj.uRef;
                % Uopt(:, k) = U;
                % % Uopt(:, k) = U


                % Uopt = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls*(obj.N_MPC));
                % Xhis(:, 2:end) = reshape(Zopt(1:obj.nStates*(obj.N_MPC+1)), obj.nStates, obj.N_MPC+1);


                U0 = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls);
                Uopt(:, k) = U0(1:obj.nControls);

                X1 = Zopt(1:(obj.N_MPC+1)*obj.nStates);
                Xhis(:, k+1) = X1(obj.nStates+1:2*obj.nStates);
                

            end

            % [Zopt, ~, ~, ~, ~, ~] = qpOASES(2*obj.G, obj.f, obj.Aeq, obj.lb, obj.ub, obj.beq, obj.beq, []);

            % Uopt = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls*(obj.N_MPC));
            % Xhis(:, 2:end) = reshape(Zopt(1:obj.nStates*(obj.N_MPC+1)), obj.nStates, obj.N_MPC+1);
        end

        % function [Xhis, Uopt] = runMPC(obj)
        %     obj.options = optimset('Display','on', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
        % 
        %     % Fix: Correct initialization
        %     Xhis = zeros(obj.nStates, obj.N_MPC+1);
        %     Uopt = zeros(obj.nControls, obj.N_MPC);
        %     Xhis(:, 1) = obj.X0;
        % 
        %     obj.beq(1:obj.nStates)
        % 
        %     % Solve the QP problem
        %     % [Zopt, ~, ~, ~, ~, ~] = qpOASES(2*obj.G, obj.f, obj.Aeq, obj.lb, obj.ub, obj.beq, obj.beq, []);
        %     [Zopt, fval, exitflag, output, lambda] = quadprog(obj.G, [], [], [], obj.Aeq, obj.beq, obj.lb, obj.ub, [], obj.options);
        % 
        %     % Fix: Ensure correct reshaping
        %     Xhis(:, 2:end) = reshape(Zopt(1:obj.nStates*(obj.N_MPC)), obj.nStates, obj.N_MPC);
        % 
        %     % Fix: Correct control input extraction
        %     Uopt = reshape(Zopt((obj.N_MPC)*obj.nStates+1:obj.nControls*obj.N_MPC+(obj.N_MPC*obj.nStates)), obj.nControls, obj.N_MPC);
        % end


        function plotResults(obj, Xhis, Uopt)
            t = 1:obj.NT+1;
            length(t)
            
            figure(1);
            subplot(3, 1, 1);
            stairs(t, Xhis(1, :), 'k', 'LineWidth', 0.5);
            grid on;
            legend('$x$', 'Interpreter', 'latex')
            title('States');
            ylabel('$x$', 'Interpreter', 'latex')
            subplot(3, 1, 2);
            stairs(t, Xhis(2, :), 'b',  'LineWidth', 0.5);
            grid on;
            legend('$z$', 'Interpreter', 'latex')
            xlabel('$Time [s]$', 'Interpreter', 'latex')
            ylabel('$z$', 'Interpreter', 'latex')
            subplot(3, 1, 3);
            stairs(t, Xhis(3, :), 'g',  'LineWidth', 0.5);
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
            
        end
    end


end

