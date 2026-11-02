% n1 n3 n5 ...
% |  |  |
% |  |  |
% |  |  |
% n2 n4 n6 ...
classdef MNA_InductorMutual < EQS_Element
    properties
        L;
        AL;
        phiDOFs;
        iDOFs;
    end
    
    methods
        function obj = MNA_InductorMutual(L)
            obj.NumberDOF = size(L, 1);
            obj.L = L;
            obj.phiDOFs = 1:2*obj.NumberDOF;
            obj.iDOFs = obj.phiDOFs(end) + (1:obj.NumberDOF);
            obj.AL = sparse([1:2:2*obj.NumberDOF, 2:2:2*obj.NumberDOF], [1:obj.NumberDOF, 1:obj.NumberDOF], ...
                    [ones(1, obj.NumberDOF), -ones(1, obj.NumberDOF)], 2*obj.NumberDOF, obj.NumberDOF);
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % L d/dti = AL'*phi;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
                Mlin(obj.iDOFs, obj.iDOFs) = obj.L;
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
            if nargout == 2
                Klin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
                Klin(obj.phiDOFs, obj.iDOFs) = obj.AL; 
                Klin(obj.iDOFs, obj.phiDOFs) = obj.AL'; 
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(3*obj.NumberDOF, 1);
            Flin = sparse(3*obj.NumberDOF, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.AL'*ys(obj.phiDOFs, :);
            obj.Solution.Current = ys(obj.iDOFs, :);
            obj.Solution.Flux = obj.L*obj.Solution.Current;
        end

        function plotSolution(obj)
            figure()
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            subplot(3,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage, LineWidth=1.5)
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            subplot(3,1,2)
            plot(obj.Solution.Time, obj.Solution.Current, LineWidth=1.5);
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
            subplot(3,1,3)
            plot(obj.Solution.Time, obj.Solution.Flux, LineWidth=1.5);
            ylabel("Flux (Wb)")
            xlabel("Time (s)")
            title("Flux")
            grid on
        end
    end
end

