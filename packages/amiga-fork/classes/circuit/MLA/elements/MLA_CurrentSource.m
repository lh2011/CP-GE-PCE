classdef MLA_CurrentSource < EQS_Element
    properties
        Current;
        B;
    end
    methods
        function obj = MLA_CurrentSource(ic)
            obj.NumberDOF = 1;
            if isnumeric(ic)
                obj.Current = @(t) ic;
            elseif isa(ic, 'function_handle')
                obj.Current = ic;
            else
                error([class(obj) ': Wrong input type!'])
            end
        end

        %%%%%%%% element equations
        % 0 = Bu
        % 0 = B'ict - i(t)
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
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
            Fnonlin(end) = -obj.Current(t);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Current = obj.B'*ys(1:obj.NumberDOFext, :);
            obj.Solution.Voltage = ys(end, :);
            obj.Solution.Power = -obj.Solution.Voltage.*obj.Solution.Current;
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

