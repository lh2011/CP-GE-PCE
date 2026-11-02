classdef MLA_Inductor < EQS_Element
    properties
        B;
        Inductance;
    end
    methods
        function obj = MLA_Inductor(L)
            obj.NumberDOF = 0;
            obj.Inductance = L;
        end

        %%%%%%%% element equations
        % B*L*B'*iLoop_prime = 0 = Bu
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            if nargout == 2
                Mlin = -obj.B*obj.Inductance*obj.B';
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Klin = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext, 1);
            Flin = sparse(obj.NumberDOFext, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Current = obj.B'*ys;
            obj.Solution.Voltage = obj.Inductance*obj.B'*ysprime;
            obj.Solution.Power = obj.Solution.Voltage.*obj.Solution.Current;
        end

        function plotSolution(obj)
            figure()
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            subplot(2,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            subplot(2,1,2)
            plot(obj.Solution.Time, obj.Solution.Current, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
        end
    end
end

