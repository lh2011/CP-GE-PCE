classdef MNA_Resistor < EQS_Element
    properties
        Resistance;
        Conductance;
    end

    methods
        function obj = MNA_Resistor(resistance, args)
            arguments
                resistance = [];
                args.CrossSection = [];
                args.Conductivity = [];
                args.Length = [];
            end
            obj.NumberDOF = 0;
            obj.calculateResistance(resistance, args.CrossSection, args.Conductivity, args.Length);
            obj.Conductance = 1/obj.Resistance;
        end

        function calculateResistance(obj, resistance, crosssection, conductivity, length)
            if ~isempty(resistance)
                obj.Resistance = resistance;
            elseif ~isempty(crosssection) && ~isempty(conductivity) && ~isempty(length)
                obj.Resistance = length/(crosssection*conductivity);
                fprintf('MNA_Resistor: Automaticalls calculated resistance to %.4g\n', obj.Resistance);
            else
                error('MNA_Resistor: Specify either resistance or crosssection, conductivity and length');
            end
        end

        %%%%%%%% element equations
        % 0 = i = AR*G*AR'*phi
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2, 2);
            Mlin = sparse(2, 2);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2, 2);
            if nargout == 2
                Klin = obj.Conductance* [1, -1; -1, 1];
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2,1);
            Flin = sparse(2,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2,2);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys;
            obj.Solution.Current = obj.Conductance*obj.Solution.Voltage;
            obj.Solution.PowerLoss = -obj.Solution.Voltage.*obj.Solution.Current;
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
            plot(obj.Solution.Time, obj.Solution.PowerLoss, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("PowerLoss (W)")
            xlabel("Time (s)")
            title("Power Loss")
            grid on
        end
    end
end

