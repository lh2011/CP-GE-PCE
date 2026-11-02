classdef MONA_Resistor < MONA_Element
    properties
        Conductivity;
        Length;
        CrossSection;

        Resistance;
        Conductance;
    end
    
    methods
        function obj = MONA_Resistor(varargin)
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
        % -A*G*A'*Psi_dot = 0
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2, 2);
            if nargout == 2
                Mlin = -obj.Conductance*[1, -1; -1, 1];
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2, 2);
            Klin = sparse(2, 2);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2, 1);
            Flin = sparse(2, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2, 2);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.A'*ysprime;
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

