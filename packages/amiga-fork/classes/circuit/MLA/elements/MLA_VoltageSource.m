classdef MLA_VoltageSource < EQS_Element
    properties
        B;
        Voltage;
    end
    methods
        function obj = MLA_VoltageSource(v)
            obj.NumberDOF = 0;
            if isnumeric(v)
                obj.Voltage = @(t) v;
            elseif isa(v, 'function_handle')
                obj.Voltage = v;
            else
                error([class(obj) ': Wrong input type!'])
            end
        end

        %%%%%%%% element equations
        % 0 = B'*u - us;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Mlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Klin = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = obj.B*obj.Voltage(t);
            Flin = sparse(obj.NumberDOFext, 1);
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Current = obj.B*ys;
            obj.Solution.Voltage = obj.Voltage(ts);
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

