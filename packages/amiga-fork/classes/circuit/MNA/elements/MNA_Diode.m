classdef MNA_Diode < EQS_Element
    properties
        Is, Vth, Rpar;
    end

    methods
        function obj = MNA_Diode(Is, Vth, Rpar)
            arguments
                Is =  1e-14;
                Vth = 2.5e-2;
                Rpar = 1e12;
            end
            obj.NumberDOF = 0;
            obj.Is = Is;
            obj.Vth = Vth;
            obj.Rpar = Rpar;
        end

        %%%%%%%% element equations
        % 0 = AD*iD
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2,2);
            Mlin = sparse(2,2);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2,2);
            Klin = sparse(2,2);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            vp = y(1);
            vm = y(2);
            I = obj.Is*(exp((vp - vm)/obj.Vth) -1) + (vp - vm)/obj.Rpar;
            Fnonlin = [I ; -I] ;
            Flin = sparse(2,1);
        end

        function J = JacobianMatrix(obj, t, y)
            vp = y(1);
            vm = y(2);
            geq = (obj.Is*exp((vp - vm)/obj.Vth))/obj.Vth + 1/obj.Rpar;
            J = [geq, -geq; -geq, geq];
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            v = [1, -1]*ys(1:2,:);
            obj.Solution.Voltage = v;
            obj.Solution.Current = obj.Is*(exp((v)/obj.Vth) -1) + (v)/obj.Rpar;
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

