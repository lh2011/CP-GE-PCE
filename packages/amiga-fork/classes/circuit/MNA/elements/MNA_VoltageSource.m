classdef MNA_VoltageSource < EQS_Element
    properties
        Voltage;
    end
    
    methods
        function obj = MNA_VoltageSource(v)
            obj.NumberDOF = 1;
            if isnumeric(v)
                obj.Voltage = @(t) ones(size(t))*v;
            elseif isa(v, 'function_handle')
                obj.Voltage = v;
            else
                error('MNA_VoltageSource: The voltage must be a scalar or function dependent on time!')
            end
        end

        %%%%%%%% element equations
        % 0 = AV*iV
        % 0 = AV'*phi - us;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(3,3);
            Mlin = sparse(3,3);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(3,3);
            if nargout == 2
                Klin = [0 0 1;
                        0 0 -1;
                        1 -1 0];
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = [0; 0; -obj.Voltage(t)];
            if nargout == 2
                Flin = sparse(3,1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(3,3);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2,:);
            obj.Solution.Current = ys(3,:);
            obj.Solution.Power = -obj.Solution.Voltage.*obj.Solution.Current;
        end

        function plotSolution(obj)
            figure()
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            subplot(3,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            subplot(3,1,2)
            plot(obj.Solution.Time, obj.Solution.Current, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
            subplot(3,1,3)
            plot(obj.Solution.Time, obj.Solution.Power, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Power (W)")
            xlabel("Time (s)")
            title("Power")
            grid on
        end
    end
end

