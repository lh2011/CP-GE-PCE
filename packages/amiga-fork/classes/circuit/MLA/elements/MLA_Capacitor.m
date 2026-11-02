classdef MLA_Capacitor < EQS_Element
    properties
        B;
        Capacitance;
    end
    methods
        function obj = MLA_Capacitor(C)
            obj.NumberDOF = 1;
            obj.Capacitance = C;
        end

        %%%%%%%% element equations
        %         0 = Bu
        % C*u_prime = B'*i_circuit
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Mlin(end, end) = obj.Capacitance;
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(1:obj.NumberDOFext, end) = obj.B;
                Klin(end, 1:obj.NumberDOFext) = obj.B';
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Current = obj.B'*ys(1:obj.NumberDOFext, :);
            obj.Solution.Voltage = ys(end, :);
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

