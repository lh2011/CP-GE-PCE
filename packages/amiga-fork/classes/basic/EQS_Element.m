classdef EQS_Element < handle
    properties (Access = public)
        Name;
        NumberDOF;
        % Store the elements, whose internal DOFS are this external DOFS
        ExternalElements = struct('Element', {}, 'LocalIndices', {});
        BoundaryConditions = struct('Boundaries', {}, 'BoundaryCondition', {});
        NumberDOFext = 0;
        Solution = struct('Time', {[]});
        G, H;           % Boundary Condition matrices
        R, K;           % Dependent and Independent degrees of freedom
        omega;          % for harmonic solutions
        % Auxiliary variables for efficiency
        LinearStiffnessMatrixShrink;
        LinearMassMatrixShrink;
        LinearForceVectorShrink;
    end

    methods (Abstract)
        % Equation system: M(t,y)y' = K(t,y)y + f(t,y)
        MassMatrix(obj, t, y)               % M: Mass matrix
        StiffnessMatrix(obj, t, y)          % K: Stiffness matrix
        ForceVector(obj, t, y)              % f: Right-hand side, which is not in K
        JacobianMatrix(obj, t, y)           % J: Nonlinear contributions of K and RHS
        postprocess(obj, ts, ys)
    end

    methods
        function addExternalElement(obj, element, indices, pos)
            if ~exist("indices", "var")
                indices = 1:element.NumberDOF;
            end
            if ~exist("pos", "var")
                pos = numel(obj.ExternalElements) + 1;
            end
            obj.ExternalElements(pos).Element = element;
            obj.ExternalElements(pos).LocalIndices = indices;
            obj.NumberDOFext = obj.NumberDOFext + numel(indices);
        end

        function resetBoundaryConditions(obj)
            obj.BoundaryConditions = struct('Boundaries', {}, 'BoundaryCondition', {});
            obj.G = []; obj.H = []; obj.R = []; obj.K = [];  % un-initialize if simulation was run before
        end

        function addBoundaryCondition(obj, boundaryCondition)
            for iBC = 1:numel(obj.BoundaryConditions)
                if boundaryCondition == obj.BoundaryConditions(iBC).BoundaryCondition
                    error("The Boundary Condition has already been added to the element!");
                end
            end
            pos = numel(obj.BoundaryConditions)+1;
            boundaryCondition.Name = inputname(2);
            obj.BoundaryConditions(pos).BoundaryCondition = boundaryCondition;
            obj.G = []; obj.H = []; obj.R = []; obj.K = [];  % un-initialize if simulation was run before
        end

        function Kred = StiffnessMatrixShrink(obj, t, yred)
            if ~exist("t", "var") || isempty(t)
                t = 0;
            end
            % update BC matrices and indices if needed
            obj.getBCmatrices();
            obj.getBCindices();
            if ~exist("yred", "var") || isempty(yred)
                yred = zeros(numel(obj.K), 1);
            end
            y = obj.reconstructSolution(t, yred);

            if isempty(obj.LinearStiffnessMatrixShrink)
                [Knonlin, Klin] = obj.StiffnessMatrix(t, y);
                obj.LinearStiffnessMatrixShrink = Klin(obj.K, obj.K) + Klin(obj.K, obj.R)*obj.G + obj.G'*Klin(obj.R, obj.K) + obj.G'*Klin(obj.R, obj.R)*obj.G;
            else
                Knonlin = obj.StiffnessMatrix(t, y);
            end
            Kred = obj.LinearStiffnessMatrixShrink + ...
                Knonlin(obj.K, obj.K) + Knonlin(obj.K, obj.R)*obj.G + obj.G'*Knonlin(obj.R, obj.K) + obj.G'*Knonlin(obj.R, obj.R)*obj.G;
        end

        function Mred = MassMatrixShrink(obj, t, yred)
            if ~exist("t", "var") || isempty(t)
                t = 0;
            end
            % update BC matrices and indices if needed
            obj.getBCmatrices();
            if ~exist("yred", "var") || isempty(yred)
                yred = zeros(numel(obj.K), 1);
            end
            y = obj.reconstructSolution(t, yred);

            if isempty(obj.LinearMassMatrixShrink)
                [Mnonlin, Mlin] = obj.MassMatrix(t, y);
                obj.LinearMassMatrixShrink = Mlin(obj.K, obj.K) + Mlin(obj.K, obj.R)*obj.G + obj.G'*Mlin(obj.R, obj.K) + obj.G'*Mlin(obj.R, obj.R)*obj.G;
            else
                Mnonlin = obj.MassMatrix(t, y);
            end
            Mred = obj.LinearMassMatrixShrink + ...
                Mnonlin(obj.K, obj.K) + Mnonlin(obj.K, obj.R)*obj.G + obj.G'*Mnonlin(obj.R, obj.K) + obj.G'*Mnonlin(obj.R, obj.R)*obj.G;
        end

        function Jred = JacobianMatrixShrink(obj, t, yred)
            if ~exist("t", "var") || isempty(t)
                t = 0;
            end
            % update BC matrices and indices if needed
            obj.getBCmatrices();
            if ~exist("yred", "var") || isempty(yred)
                yred = zeros(numel(obj.K), 1);
            end
            y = obj.reconstructSolution(t, yred);
            J = obj.JacobianMatrix(t, y);

            % JacobianMatrix(t, y) does not consider StiffnessMatrix from
            % product rule, Shrink does
            Jred = obj.StiffnessMatrixShrink(t, yred) + ...
                J(obj.K, obj.K) + J(obj.K, obj.R)*obj.G + obj.G'*J(obj.R, obj.K) + obj.G'*J(obj.R, obj.R)*obj.G;
        end
        
        function Fred = ForceVectorShrink(obj, t, yred)
            if ~exist("t", "var") || isempty(t)
                t = 0;
            end
            % update BC matrices and indices if needed
            obj.getBCmatrices();
            if ~exist("yred", "var") || isempty(yred)
                yred = zeros(numel(obj.K), 1);
            end
            y = obj.reconstructSolution(t, yred);

            if isempty(obj.LinearForceVectorShrink)
                [Fnonlin, Flin] = obj.ForceVector(t, y);
                obj.LinearForceVectorShrink = Flin(obj.K) + obj.G'*Flin(obj.R);
            else
                Fnonlin = obj.ForceVector(t, y);
            end

            Fred = obj.LinearForceVectorShrink + Fnonlin(obj.K) + obj.G'*Fnonlin(obj.R);

            % additionaly considerations of nonzero or time varying boundary conditions
            [b, db] = obj.getBCvalues(t);
            if any(b~=0, "all")
                [Knonlin, Klin] = obj.StiffnessMatrix(t, y);
                Kmat = Knonlin + Klin;
                Fred = Fred + (Kmat(obj.K, obj.R) + obj.G'*Kmat(obj.R, obj.R))*obj.H*b;
            end
            if any(db~=0, "all")
                [Mnonlin, Mlin] = obj.MassMatrix(t, y);
                Mmat = Mnonlin + Mlin;
                Fred = Fred - obj.G'*Mmat(obj.R,obj.R)*obj.H*db;
            end
        end
        
        function [ys, ysprime] = reconstructSolution(obj, ts, yred, yredprime)
            ys = zeros(obj.NumberDOF, numel(ts));
            bs = zeros(numel(obj.R), numel(ts));
            if nargout == 1
                for ti = 1:numel(ts)
                    bs(:, ti) = obj.getBCvalues(ts(ti));
                end
                ys(obj.K, :) = yred;
                ys(obj.R, :) = obj.G*yred + obj.H*bs;
            else
            % reconstruct also derivatives
                ysprime = zeros(obj.NumberDOF, numel(ts));
                bsprime = zeros(numel(obj.R), numel(ts));
                for ti = 1:numel(ts)
                    [bs(:, ti), bsprime(:, ti)] = obj.getBCvalues(ts(ti));
                end
                ys(obj.K, :) = yred;
                ys(obj.R, :) = obj.G*yred + obj.H*bs;
                ysprime(obj.K, :) = yredprime;
                ysprime(obj.R, :) = obj.G*yredprime + obj.H*bsprime;
            end
        end

        function ysred = shrinkSolution(obj, ys)
            ysred = ys(obj.K, :);
        end

        function setOmega(obj, omega)
            obj.omega = omega;
        end

        function W = OmegaMatrix(obj)
            W = obj.omega * speye(obj.NumberDOF);
        end

        function [G, H] = getBCmatrices(obj)
            % values were already called and initialized
            if issparse(obj.G) && issparse(obj.H)
                G = obj.G;
                H = obj.H;
                return
            end
            G = sparse(obj.NumberDOF, obj.NumberDOF);       % G = (RxK)
            H = sparse(obj.NumberDOF, obj.NumberDOF);       % H = (RxR)
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                [G_iel, H_iel] = BC.getBCmatrices();
                [localDofsR, localDofsK] = BC.getBCindices();
                G(localDofsR, localDofsK) = G(localDofsR, localDofsK) + G_iel;
                H(localDofsR, localDofsR) = H(localDofsR, localDofsR) + H_iel;
            end
            [dofsR, dofsK] = obj.getBCindices();
            G = G(dofsR, dofsK);
            H = H(dofsR, dofsR);
            obj.G = G;
            obj.H = H;
        end

        function [dofsR, dofsK] = getBCindices(obj)
            % values were already called and initialized
            if ~isempty(obj.K) || ~isempty(obj.R)
                dofsR = obj.R;
                dofsK = obj.K;
                return
            end
            % empty BC, so all indices will be kept
            if isempty(obj.BoundaryConditions)
                dofsR = zeros(0,1);
                dofsK = reshape(1:obj.NumberDOF, obj.NumberDOF, 1);
                obj.R = dofsR;
                obj.K = dofsK;
                return
            end
            % otherwise: collect all indices which will be removed
            dofsR = [];
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                localDofsR = BC.getBCindices();
                dofsR = reshape(union(dofsR, localDofsR), 1, []);
            end
            dofsK = setdiff((1:obj.NumberDOF), dofsR);

            obj.R = dofsR;
            obj.K = dofsK;
            assert(numel(dofsR) == numel(unique(dofsR)))
            assert(numel(dofsK) == numel(unique(dofsK)))
            assert(numel(dofsK)+numel(dofsR) == obj.NumberDOF);
        end

        function [b, db] = getBCvalues(obj, t)
            b = sparse(obj.NumberDOF, 1);
            db = sparse(obj.NumberDOF, 1);
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                localDofsR = BC.getBCindices();
                [b_iel, db_iel] = BC.getBCvalues(t);
                b(localDofsR) = b(localDofsR) + b_iel;
                db(localDofsR) = db(localDofsR) + db_iel;
            end
            dofsR = obj.getBCindices();
            b = b(dofsR, 1);
            db = db(dofsR, 1);
        end

    end
end
