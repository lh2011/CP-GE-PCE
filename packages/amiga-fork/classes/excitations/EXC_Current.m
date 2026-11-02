classdef EXC_Current < EXC
    properties (Access = public)
        CircuitInterface;
    end

    methods
        function obj = EXC_Current(material, phaseDefinition, currentFunction, args)
            arguments
                material MAT
                phaseDefinition 
                currentFunction 
                args.NumberWindings = 1;
                args.PlotColor = [0,0,0];
            end
            obj@EXC(material, horzcat(phaseDefinition.Regions));
            material.addCurrentPostprocessing(obj, obj.RegionNumbers);

            obj.ExcitationFunction = currentFunction;
            obj.CircuitInterface = MNA_CircuitInterface(material, phaseDefinition, "NumberWindings", args.NumberWindings, "PlotColor", args.PlotColor);
            obj.IsSpaceDependent = false;
            if isa(currentFunction, "function_handle")
                obj.IsTimeDependent = true;
            else
                obj.IsTimeDependent = false;
            end
            obj.updateExcitation()
        end

        function updateExcitation(obj)
            X = obj.CircuitInterface.getX();
            % reduction for active indices if excitation is not set for all regions which belong to the material
            [~, ~, indices] = obj.ParentMaterial.extractMeshSpace(horzcat(obj.CircuitInterface.Phases.Regions));
            obj.RHS = X(indices);
        end

        function J = getCurrentDensity(obj, t, x, ~, ~, reg)
            if obj.IsTimeDependent
                current = obj.ExcitationFunction(t);
            else
                current = obj.ExcitationFunction;
            end
            % find correct phase
            for iPhase = 1:numel(obj.CircuitInterface.Phases)
                nWind = obj.CircuitInterface.NumberWindings;
                if ismember(reg, obj.CircuitInterface.Phases(iPhase).Regions)
                    area = obj.CircuitInterface.Phases(iPhase).Area;
                    direction = obj.CircuitInterface.Phases(iPhase).Direction;
                    J = direction*nWind*current/area * ones(size(x));
                    break;
                end
            end
        end

        function Pv = getLossDensity(obj, t, x, ~, ~, reg)
            Pv = zeros(size(x));
        end
    end
end