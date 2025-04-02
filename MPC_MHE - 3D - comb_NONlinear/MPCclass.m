


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
        zPrev
        dt
        lb
        ub
        xRef
        lbuRef
        ubuRef
        options
    end

    methods
        function obj = MPCclass(N_MPC, Ac, Bc, X0, dt, lb, ub, Q, R, nStates, nControls,options, xRef, lbuRef, ubuRef)
            obj.N_MPC = N_MPC;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.X0 = X0;
            obj.dt = dt;
            obj.nStates = nStates;
            obj.nControls = nControls;
            obj.Q = Q;
            obj.R = R;
            obj.lb = lb;
            obj.ub = ub;
            obj.xRef = xRef;
            obj.lbuRef = lbuRef;
            obj.ubuRef = ubuRef;
            obj.options = options;
            obj = obj.initialize();
            
        end

        % Function for discretezise A and B matrices
        function obj = initialize(obj)
            obj.A = expm(obj.Ac * obj.dt);
            obj.B = (expm(obj.Ac * obj.dt) - eye(size(obj.Ac))) * (obj.Ac \ obj.Bc);
            
            [~,P,~] = dlqr(obj.A, obj.B, obj.Q, obj.R);
            obj.QN = P;
            obj.zPrev = [obj.X0;zeros(obj.N_MPC*obj.nStates,1);zeros((obj.N_MPC)*obj.nControls,1)];

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

        function uRef = computeReferenceInput(obj)

            GU = 10*eye(obj.nControls);
            fU = zeros(obj.nControls, 1);

            AeqU = obj.B;
            beqU = (eye(obj.nStates) - obj.A) * obj.xRef;

            options = optimset('Display','on', 'Diagnostics','on', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
            
            uRef = quadprog(GU, fU, [], [], AeqU, beqU, obj.lbuRef, obj.ubuRef, [], options);
            
        end

        function [xNextTrue, Uopt] = runMPC(obj,xCurrent)
        
        % uRef = pinv(obj.B) * ((eye(obj.nStates) - obj.A) * obj.xRef);
        obj.beq(1:obj.nStates) = xCurrent;
        [Zopt, ~, ~, ~, ~] = quadprog(2*obj.G, obj.f, [],[],obj.Aeq, obj.beq,obj.lb,obj.ub, obj.zPrev, obj.options);
        obj.zPrev=Zopt;

        U0 = Zopt((obj.N_MPC+1)*obj.nStates+1 : (obj.N_MPC+1)*obj.nStates + obj.nControls);
        Uopt = U0(1:obj.nControls);

        X1 = Zopt(1:(obj.N_MPC+1)*obj.nStates);
        xNextTrue=X1(obj.nStates+1:2*obj.nStates);


        end



        function plotResults(obj, Xopt, Uopt, xRef)

            % Plotting the results
            %t = 0:obj.NT;
            %length(t)
            
            figure(1);
            subplot(3, 1, 1);
            stairs(Xopt(1, :), 'k', 'LineWidth', 1);hold on
            plot([0, size(Xopt,2)],[xRef(1), xRef(1)])
            %text(30, 0.9*xRef(1), 'xRef', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
            grid on;
            ylim([1.1*min(Xopt(1,:))-1e-4, 1.1*max(Xopt(1,:))+1e-4])
            legend('$x$', 'Interpreter', 'latex')
            title('States');
            ylabel('$x$', 'Interpreter', 'latex')

            subplot(3, 1, 2);
            stairs(Xopt(2, :), 'b',  'LineWidth', 1);hold on
            plot([0, size(Xopt,2)],[xRef(2), xRef(2)])
            %text(30, 0.9*xRef(2), 'zRef', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
            grid on;
            ylim([1.1*min(Xopt(2,:))-1e-4, 1.1*max(Xopt(2,:))+1e-4])
            legend('$z$', 'Interpreter', 'latex')
            ylabel('$z$', 'Interpreter', 'latex')

            subplot(3, 1, 3);
            stairs(Xopt(3, :),'color', [0.9290 0.6940 0.1250],  'LineWidth', 1);hold on
            plot([0, size(Xopt,2)],[xRef(3), xRef(3)])
            %text(30, 0.9*xRef(3), 'thetaRef', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
            grid on;
            ylim([1.1*min(Xopt(3,:))-1e-4, 1.1*max(Xopt(3,:))+1e-4])
            legend('$\theta$', 'Interpreter', 'latex')
            xlabel('$Iterations$', 'Interpreter', 'latex')
            ylabel('$\theta$', 'Interpreter', 'latex')


            figure(2);
            subplot(2, 1, 1)
            stairs(Uopt(1,:), 'r.-',  'LineWidth', 0.5);
            grid on;
            legend('$ux$', 'Interpreter', 'latex');
            ylabel('$ux$', 'Interpreter', 'latex');
            title("Control inputs")
            subplot(2, 1, 2)
            stairs(Uopt(2,:), 'r.-',  'LineWidth', 0.5);
            grid on;
            legend('$uz$', 'Interpreter', 'latex');
            ylabel('$uz$', 'Interpreter', 'latex');
            xlabel('Iterations','Interpreter','latex');


            %figure(3)
            %plot(Xopt(1, :), Xopt(2, :));
            %grid on;
            %ylabel('$x$', 'Interpreter', 'latex');
            %xlabel('$z$',  'Interpreter', 'latex');
            %title('PhasePortrait of x and z')
        end
    end


end

