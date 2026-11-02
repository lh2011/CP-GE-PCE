classdef MNA_StrandedConductor < EQS_Element & MNA_CircuitInterface
    properties
        Circuit;
        aDOFs, phiDOFs, iDOF;
        X;
        Rcoil;
    end

    methods
        function obj = MNA_StrandedConductor(material, phaseDefinition, circuit, nodes, args)
            arguments
                material MAT
                phaseDefinition 
                circuit MNA
                nodes
                args.NumberWindings = 1;
                args.PlotColor = [0, 0, 0];
                args.Resistance = [];
                args.WireDiameter = [];
                args.FillFactor = [];
                args.Conductivity = [];
            end
            obj@MNA_CircuitInterface(material, phaseDefinition, "NumberWindings", args.NumberWindings, "PlotColor", args.PlotColor);
            material.addCurrentPostprocessing(obj, horzcat(phaseDefinition.Regions));
            obj.calculateResistance(args.Resistance, args.WireDiameter, args.FillFactor, args.Conductivity);
            if obj.Sigma ~= 0
                material.setConductivity(0);
                obj.Sigma = 0;
                warning('MNA_StrandedConductor: The conductivity of the connected material is set to 0');
            end
            obj.NumberDOF = 1;
            obj.addExternalElement(material.Region, material.Space.gnum);
            obj.Circuit = circuit;
            obj.X = obj.getX();
            
            nodeIndices = circuit.Nodes.getNodeIndices(nodes);
            obj.addExternalElement(circuit.Nodes, nodeIndices);
            % Dofs ordering
            obj.aDOFs = 1:material.Space.ndof;
            obj.phiDOFs = material.Space.ndof + (1:2);
            obj.iDOF = material.Space.ndof + 3;
        end

        % Field:          (M_sigma*a')  = (-K_nu*a)   + A*i
        % Circuit:        (M_crc*phi')  = (K_crc*phi) + A*i
        % Coupling:              X.T*a' = A.T*phi - R*i
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Mlin(obj.iDOF, obj.aDOFs) = obj.X';
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.aDOFs, obj.iDOF) = obj.X;
                Klin(obj.phiDOFs, obj.iDOF) = obj.A;
                Klin(obj.iDOF, obj.phiDOFs) = obj.A';
                Klin(obj.iDOF, obj.iDOF) = -obj.Rcoil;
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
        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.A'*ys(obj.phiDOFs, :);
            obj.Solution.Current = ys(obj.iDOF, :);
            obj.Solution.FluxLinkage = obj.X'*ys(obj.aDOFs, :);
            obj.Solution.PowerLoss = -obj.Rcoil*obj.Solution.Current.^2;
        end

        function J = getCurrentDensity(obj, t, x, ~, ~, reg)
            if isscalar(obj.Solution.Time)
                current = obj.Solution.Current;
            else
                current = interp1(obj.Solution.Time, obj.Solution.Current, t);
            end
            % find correct phase
            for iPhase = 1:numel(obj.Phases)
                if ismember(reg, obj.Phases(iPhase).Regions)
                    area = obj.Phases(iPhase).Area;
                    direction = obj.Phases(iPhase).Direction;
                    J = direction*obj.NumberWindings*current/area * ones(size(x));
                    break;
                end
            end
        end

        function calculateResistance(obj, resistance, wireDiameter, fillFactor, sigma)
            if ~isempty(resistance)
                obj.Rcoil = resistance;
                fprintf('MNA_StrandedConductor: Manually set phase resistance to %.4g\n', obj.Rcoil);
            elseif ~isempty(wireDiameter) && ~isempty(sigma)
                wireCrossSection = wireDiameter^2*pi/4;
                nPhases = numel(obj.Phases);
                obj.Rcoil = obj.ParentRegion.Length*nPhases*obj.NumberWindings/(sigma*wireCrossSection);
                fprintf(['MNA_StrandedConductor: Automatically calculated phase resistance with wire ' ...
                    'diameter %.4g, and conductivity %.4g to be %.4g\n'], wireDiameter, sigma, obj.Rcoil);
            elseif ~isempty(fillFactor) && ~isempty(sigma)
                phaseArea = obj.Phases(1).Area;
                wireCrossSection = phaseArea*fillFactor/obj.NumberWindings;
                nPhases = numel(obj.Phases);
                obj.Rcoil = obj.ParentRegion.Length*nPhases*obj.NumberWindings/(sigma*wireCrossSection);
                fprintf(['MNA_StrandedConductor: Automatically calculated phase resistance with fill ' ...
                    'factor %.4g, wire diameter %.4g, and conductivity %.4g to be %.4g\n'], fillFactor, sqrt(wireCrossSection*4/pi), sigma, obj.Rcoil);
            else
                obj.Rcoil = 0;
                fprintf('MNA_StrandedConductor: No further resistance specified. Automatically set phase resistance to %.4g\n', obj.Rcoil);
            end
        end

        function Pv = getLossDensity(obj, t, x, ~, ~, reg)
            if isscalar(obj.Solution.Time)
                current = obj.Solution.Current;
            else
                current = interp1(obj.Solution.Time, obj.Solution.Current, t);
            end
            % find correct phase
            for iPhase = 1:numel(obj.Phases)
                if ismember(reg, obj.Phases(iPhase).Regions)
                    area = obj.Phases(iPhase).Area;
                    ReqRegion = obj.Rcoil/numel(obj.Phases); % equivalent resistance for part of phase
                    Pv = ReqRegion*current.^2/(area*obj.ParentRegion.Length) * ones(size(x));
                    break;
                end
            end
        end

        function plotSolution(obj)
            figure()
            sgtitle(obj.Name)
            subplot(3,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage, "Color", TUDa_getColor("1b"), LineWidth=1.5)
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
            plot(obj.Solution.Time, obj.Solution.FluxLinkage, "Color", TUDa_getColor("1b"), LineWidth=1.5)
            ylabel("Magnetic Flux Linkage (Wb)")
            xlabel("Time (s)")
            title("Flux Linkage")
            grid on
        end
    end
end


