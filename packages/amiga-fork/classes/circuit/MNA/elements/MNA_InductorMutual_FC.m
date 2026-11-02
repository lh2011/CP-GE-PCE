% Flux charge implementation for a general (linear) inductor with n inputs
% n1 n3 n5 ...
% |  |  |
% |  |  |
% |  |  |
% n2 n4 n6 ...
classdef MNA_InductorMutual_FC < EQS_Element
    properties
        L;
        AL;
        NumberPorts;
        phiDOFs;    % Dofs for the voltage potentials
        PhiDOFs;    % Dofs for the magnetic fluxes
        iDOFs;      % Dofs for the currents
    end
    methods
        function obj = MNA_InductorMutual_FC(L)
            obj.NumberPorts = size(L, 1);
            obj.NumberDOF = 2*size(L, 1);
            obj.phiDOFs = 1:obj.NumberDOF;
            obj.PhiDOFs = obj.phiDOFs(end) + (1:obj.NumberPorts);
            obj.iDOFs = obj.PhiDOFs(end) + (1:obj.NumberPorts);
            obj.L = L;
            obj.AL = sparse([1:2:obj.NumberDOF, 2:2:obj.NumberDOF], [1:obj.NumberPorts, 1:obj.NumberPorts], ...
                    [ones(1, obj.NumberPorts), -ones(1, obj.NumberPorts)], obj.NumberDOF, obj.NumberPorts);
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % d/dt Phi = AL'*phi = uL
        % 0 = Phi - Li;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2*obj.NumberDOF, 2*obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(2*obj.NumberDOF, 2*obj.NumberDOF);
                Mlin(obj.PhiDOFs, obj.PhiDOFs) = eye(obj.NumberPorts);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2*obj.NumberDOF, 2*obj.NumberDOF);
            if nargout == 2
                Klin = sparse(2*obj.NumberDOF, 2*obj.NumberDOF);
                Klin(obj.iDOFs, obj.iDOFs) = -obj.L;
                Klin(obj.iDOFs, obj.PhiDOFs) = eye(obj.NumberPorts);
                Klin(obj.phiDOFs, obj.iDOFs) = obj.AL; 
                Klin(obj.PhiDOFs, obj.phiDOFs) = obj.AL'; 
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2*obj.NumberDOF, 1);
            Flin = sparse(2*obj.NumberDOF, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2*obj.NumberDOF, 2*obj.NumberDOF);
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.AL'*ys(obj.phiDOFs, :);
            obj.Solution.Flux = ys(obj.PhiDOFs, :);
            obj.Solution.Current = ys(obj.iDOFs, :);
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

