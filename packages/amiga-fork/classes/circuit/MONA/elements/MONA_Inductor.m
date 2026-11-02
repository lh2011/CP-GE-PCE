classdef MONA_Inductor < MONA_Element
    properties
        Inductance;
    end
    
    methods
        function obj = MONA_Inductor(L)
            obj.NumberDOF = 0;
            obj.Inductance = L;
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % 0 = AL*L^-1AL'*Psi
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2,2);
            Mlin = sparse(2,2);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2,2);
            Klin = obj.Inductance^-1* [1, -1; -1, 1];
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
            obj.Solution.Voltage = obj.A'*ysprime;
            obj.Solution.Current = obj.Inductance^-1 *obj.A'*ys;
            obj.Solution.Power = obj.Solution.Voltage.*obj.Solution.Current;
            obj.Solution.Flux = obj.A'*ys;
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

