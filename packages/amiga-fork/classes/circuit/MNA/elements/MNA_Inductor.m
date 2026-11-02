classdef MNA_Inductor < EQS_Element
    properties
        Inductance;
    end
    
    methods
        function obj = MNA_Inductor(L)
            obj.NumberDOF = 1;
            obj.Inductance = L;
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % L d/dt i = AL'*phi = uL
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(3,3);
            Mlin = [0 0 0; 
	                0 0 0; 
	                0 0 obj.Inductance];
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(3,3);
            Klin = [0 0 1; 
                    0 0 -1; 
                    1 -1 0];
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(3,1);
            Flin = sparse(3,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(3,3);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2, :);
            obj.Solution.Current = ys(3, :);
            obj.Solution.Power = obj.Solution.Voltage.*obj.Solution.Current;
            obj.Solution.Flux = obj.Inductance*obj.Solution.Current;
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
            plot(obj.Solution.Time, obj.Solution.Flux, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Flux (Wb)")
            xlabel("Time (s)")
            title("Flux")
            grid on
        end
    end
end

