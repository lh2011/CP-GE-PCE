classdef MONA_SolidConductor < EQS_Element & MNA_CircuitInterface
    properties
        Circuit;
        aDOFs, psiDOFs, qDOF;
        Conductance;
        Q;
        Msigma;
    end

    methods
        function obj = MONA_SolidConductor(material, phaseDefinition, circuit, nodes, args)
            arguments
                material MAT
                phaseDefinition 
                circuit MONA
                nodes 
                args.PlotColor = [0, 0, 0];
            end
            obj@MNA_CircuitInterface(material, phaseDefinition, "PlotColor", args.PlotColor);
            material.addCurrentPostprocessing(obj, horzcat(phaseDefinition.Regions));

            if obj.Sigma == 0
                warning("The conductivity of the solid conductor material must be > 0!");
            end
            obj.NumberDOF = 1;
            
            obj.addExternalElement(material.Region, material.Space.gnum);
            obj.Circuit = circuit;
            obj.Q = obj.getQ();
            obj.Msigma = obj.getMsigma();
            
            nodeIndices = circuit.Nodes.getNodeIndices(nodes);
            obj.addExternalElement(circuit.Nodes, nodeIndices);

            obj.aDOFs = 1:material.Space.ndof;
            obj.psiDOFs = material.Space.ndof + (1:2);
            obj.qDOF = material.Space.ndof + 3;
            obj.Conductance = sum(vertcat(obj.Phases.Conductance));
        end

        % Field:        M_sigma*a' - Q*A.T*psi'  = (-K_nu*a)
        % Circuit:      (M_crc*phi') - A*q'  = (K_crc*phi)
        % Coupling:     Q.T*a' + q' - G*A.T*psi'  = 0
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Mlin(obj.aDOFs, obj.psiDOFs) = - obj.Q*obj.A';
                Mlin(obj.psiDOFs, obj.qDOF)  = - obj.A;
                Mlin(obj.qDOF, obj.aDOFs)    = obj.Q';
                Mlin(obj.qDOF, obj.psiDOFs)  = - obj.Conductance*obj.A';
                Mlin(obj.qDOF, obj.qDOF)     = 1;
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.A'*ysprime(obj.psiDOFs, :);
            obj.Solution.Current = ysprime(obj.qDOF, :);
            obj.Solution.PowerLoss = zeros(1, numel(ts));
            obj.Solution.AzPrime = ysprime(obj.aDOFs, :);
            for iTime = 1:numel(ts)
                usol = obj.Solution.Voltage(iTime);
                adot = ysprime(obj.aDOFs, iTime);
                obj.Solution.PowerLoss(iTime) = -adot'*obj.Msigma*adot + 2*usol*obj.Q'*adot - usol*obj.Conductance*usol;
            end
        end

        function J = getCurrentDensity(obj, t, x, ~, AzPrime, reg)
            if isscalar(obj.Solution.Time)
                deltaPhi = obj.Solution.Voltage / obj.ParentRegion.Length;
            else
                deltaPhi = interp1(obj.Solution.Time, obj.Solution.Voltage, t) / obj.ParentRegion.Length;
            end
            % find correct phase
            for iPhase = 1:numel(obj.Phases)
                if ismember(reg, obj.Phases(iPhase).Regions)
                    direction = obj.Phases(iPhase).Direction;
                    J = direction*obj.Sigma*(-AzPrime + deltaPhi);
                    break;
                end
            end
        end

        function Pv = getLossDensity(obj, t, x, ~, AzPrime, reg)
            if isscalar(obj.Solution.Time)
                deltaPhi = obj.Solution.Voltage / obj.ParentRegion.Length;
            else
                deltaPhi = interp1(obj.Solution.Time, obj.Solution.Voltage, t) / obj.ParentRegion.Length;
            end
            % find correct phase
            for iPhase = 1:numel(obj.Phases)
                if ismember(reg, obj.Phases(iPhase).Regions)
                    direction = obj.Phases(iPhase).Direction;
                    Pv = direction*obj.Sigma*(-AzPrime + deltaPhi).^2;
                    break;
                end
            end
        end

        function plotSolution(obj)
            figure()
            sgtitle(obj.Name)
            subplot(2,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage, "Color", TUDa_getColor("1b"), LineWidth=1.5)
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            subplot(2,1,2)
            plot(obj.Solution.Time, obj.Solution.Current, "Color", TUDa_getColor("1b"), LineWidth=1.5)
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
        end
    end
end


