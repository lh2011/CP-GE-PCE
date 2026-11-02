classdef MNA_Capacitor < EQS_Element
    properties
        Capacitance;
    end
    
    methods
        function obj = MNA_Capacitor(C)
            obj.NumberDOF = 0;
            obj.Capacitance = C;
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % d/dt Phi = AL'*phi = uL
        % 0 = Phi - Li;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2,2);
            Mlin = [-obj.Capacitance ,obj.Capacitance; 
                    obj.Capacitance, -obj.Capacitance];
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2,2);
            Klin = sparse(2,2);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2,1);
            Flin = sparse(2,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2,2);
        end

        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2, :);
            obj.Solution.Current = obj.Capacitance* [1, -1]*ysprime(1:2, :);
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

