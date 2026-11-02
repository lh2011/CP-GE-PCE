classdef MLA_Loops < EQS_Element

    methods
        function obj = MLA_Loops()
            obj.NumberDOFext = 0;
            obj.NumberDOF = 0;
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
            obj.Solution.LoopCurrents = ys;
        end
        
        function plotSolution(obj)
            figure();
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            plot(obj.Solution.Time, obj.Solution.LoopCurrents)
            title("Loop Currents")
            ylabel("Current (A)")
            xlabel("Time (s)");
            legend("Loop "+string(obj.LoopNumbers));
            grid on
        end
    end
end

