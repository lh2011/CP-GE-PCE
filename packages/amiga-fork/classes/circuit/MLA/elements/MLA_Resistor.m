classdef MLA_Resistor < EQS_Element
    properties
        B;
        
        Conductivity;
        Length;
        CrossSection;

        Resistance;
        Conductance;
    end
    methods
        function obj = MLA_Resistor(varargin)
            obj.NumberDOF = 0;
            if nargin == 1
                obj.Resistance = varargin{1};
            elseif nargin == 3
                obj.Conductivity = varargin{1};
                obj.Length = varargin{2};
                obj.CrossSection = varargin{3};
                obj.Resistance = obj.CrossSection/(obj.Conductivity*obj.Length);
            else
                error([class(obj) ': Wrong number of input arguments!'])
            end
            obj.Conductance = 1/obj.Resistance;
        end

        %%%%%%%% element equations
        % 0 = Bu = B*R*B'*iLoop
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Mlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            if nargout == 2
                Klin = obj.B*obj.Resistance*obj.B';
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext, 1);
            Flin = sparse(obj.NumberDOFext, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Current = obj.B'*ys;
            obj.Solution.Voltage = obj.Resistance*obj.Solution.Current;
            obj.Solution.PowerLoss = -obj.Solution.Voltage.*obj.Solution.Current;
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

