classdef MNA_SolidConductor < EQS_Element & MNA_CircuitInterface
    properties
        Circuit;
        aDOFs, phiDOFs, iDOF;
        Conductance;
        Q;
        Msigma;
    end

    methods
        function obj = MNA_SolidConductor(material, phaseDefinition, circuit, nodes, args)
            arguments
                material MAT
                phaseDefinition 
                circuit MNA
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
            obj.phiDOFs = material.Space.ndof + (1:2);
            obj.iDOF = material.Space.ndof + 3;
            obj.Conductance = sum(vertcat(obj.Phases.Conductance));
        end

        % Field:        M_sigma*a'   = (-K_nu*a) + Q*A.T*phi
        % Circuit:      (M_crc*phi'  = (K_crc*phi) + A*i 
        % Coupling:     Q.T*a'       = G*A.T*phi - i
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Mlin(obj.iDOF, obj.aDOFs) = obj.Q';
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.aDOFs, obj.phiDOFs) = obj.Q*obj.A';
                Klin(obj.phiDOFs, obj.iDOF)  = obj.A;
                Klin(obj.iDOF, obj.phiDOFs)  = obj.Conductance*obj.A';
                Klin(obj.iDOF, obj.iDOF)     = -1;
            end
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
            obj.Solution.Voltage = obj.A'*ys(obj.phiDOFs, :);
            obj.Solution.Current = ys(obj.iDOF, :);
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

        % function checkCurrentDensity(obj)
        %     % I = int (sigma(-Adot + deltaPhi)) dA
        %     Azdot = obj.Solution.AzPrime;
        %     [submesh, subspace, subindices] = obj.ParentMaterial.extractMeshSpace(obj.Phases.Regions);
        %     Azdot = Azdot(subindices, :);
        %     azel = zeros([size(subspace.connectivity), size(Azdot, 2)]);
        %     for iel = 1:submesh.nel
        %         azel(:, iel, :) = Azdot(subspace.connectivity(:, iel), :);
        %     end
        % 
        %     deltaPhi = obj.Solution.Voltage/obj.ParentRegion.Length;
        %     quad = submesh.jacdet.*submesh.quad_weights;
        %     NiAi = sum(subspace.shape_functions.* permute(azel, [4,1,2,3]), 2);
        %     I1 = permute(obj.Sigma*sum(-permute(quad, [1,3,2]).*NiAi, [1,2,3]), [1,4,2,3]);
        %     I2 = obj.Sigma*sum(quad, "all").* deltaPhi;
        %     % I =  + permute(deltaPhi, [1,3,4,2])), [2,3]);
        %     % I = permute(I, [1,4,2,3]);
        %     Q = zeros(obj.ParentMaterial.Space.ndof, 1);
        %     Q(subindices) = obj.Sigma*op_f_v1(subspace, submesh);
        % 
        %     II1 = -obj.Q' * obj.Solution.AzPrime;
        %     II2 = obj.Conductance*obj.Solution.Voltage;
        %     I = I1 + I2;
        % 
        %     % figure
        %     % hold on
        %     % plot(obj.Solution.Time, obj.Solution.Current);
        %     % plot(obj.Solution.Time, I);
        %     figure
        %     hold on
        %     plot(obj.Solution.Time, I1+I2);
        %     plot(obj.Solution.Time, II1+II2);
        %     plot(obj.Solution.Time, obj.Solution.Current);
        % end

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


