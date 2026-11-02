classdef MNA_Nodes < EQS_Element
    properties
        NodeNumbers = [];
    end

    methods
        function obj = MNA_Nodes()
            obj.NumberDOFext = 0;
            obj.NumberDOF = 0;
            BC_Dirichlet(obj, 1); % zeroPotential, automatically added to node with first index
        end

        function [indNodes, indNewNodes] = addNodes(obj, nodes)
            newNodes = setdiff(nodes, obj.NodeNumbers, 'stable');
            obj.NodeNumbers = [obj.NodeNumbers, newNodes];
            obj.NumberDOF = numel(obj.NodeNumbers);
            indNodes = obj.getNodeIndices(nodes);
            indNewNodes = obj.getNodeIndices(newNodes);
            % Update Zero-Potential Condition for new independent DOFs
            obj.BoundaryConditions(1).BoundaryCondition.KeepDOFs = setdiff(1:obj.NumberDOF, obj.BoundaryConditions(1).BoundaryCondition.RemoveDOFs);
            obj.BoundaryConditions(1).BoundaryCondition.CK = sparse(numel(obj.BoundaryConditions(1).BoundaryCondition.RemoveDOFs), numel(obj.BoundaryConditions(1).BoundaryCondition.KeepDOFs));
        end

        function setZeroPotentialNode(obj, nodeNumber)
            rDofs = obj.getNodeIndices(nodeNumber);
            kDofs = setdiff(1:obj.NumberDOF, rDofs);
            obj.BoundaryConditions(1).BoundaryCondition.RemoveDOFs = rDofs;
            obj.BoundaryConditions(1).BoundaryCondition.KeepDOFs = kDofs;
        end

        function nodesIndices = getNodeIndices(obj, nodeNumber)
            [~, ~, nodesIndices] = intersect(nodeNumber, obj.NodeNumbers, 'stable');
            nodesIndices = reshape(nodesIndices, 1, []);
        end

        %%%%%%%% element equations: Nodes have no contribution
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
            Mlin = sparse(obj.NumberDOF, obj.NumberDOF);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOF, obj.NumberDOF);
            Klin = sparse(obj.NumberDOF, obj.NumberDOF);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOF, 1);
            Flin = sparse(obj.NumberDOF, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOF, obj.NumberDOF);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.VoltagePotenials = ys;
        end
        
        function plotSolution(obj)
            figure();
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            plot(obj.Solution.Time, obj.Solution.VoltagePotenials)
            title("Nodes Voltage Potential")
            ylabel("Voltage (V)")
            xlabel("Time (s)");
            legend("Node "+string(obj.NodeNumbers));
            grid on
        end
    end
end

