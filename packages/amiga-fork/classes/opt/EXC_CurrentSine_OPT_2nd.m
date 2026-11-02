classdef EXC_CurrentSine_OPT_2nd < EXC_CurrentSine_OPT & OPT_Element_2nd
    properties (Access = public)
        d2Current_dIAppdPhase0;
        d2Current_dPhase02;
    end

    methods
        function obj = EXC_CurrentSine_OPT_2nd(material, phasesStruct, iPeak, omega, uvwshift, phase0, args)
            arguments
                material 
                phasesStruct 
                iPeak 
                omega 
                uvwshift 
                phase0
                args.NumberWindings = 1;
                args.PlotColor = [0,0,0];
            end
           obj@EXC_CurrentSine_OPT(material, phasesStruct, iPeak, omega, uvwshift, phase0, NumberWindings=args.NumberWindings, PlotColor=args.PlotColor);
           obj.d2Current_dIAppdPhase0 = @(t) cos(obj.Omega*t + obj.UVWshift + obj.Phase0);
           obj.d2Current_dPhase02 = @(t) -obj.IApp*sin(obj.Omega*t + obj.UVWshift + obj.Phase0);
        end
        
        function Lxx = function_Lxx(obj, t, x, y, z)
            Lxx = zeros(obj.ParentRegion.NumberOptimizationParameters, obj.ParentRegion.NumberOptimizationParameters, obj.ParentRegion.NumberOptimizationResponses);
            if obj.ParentRegion.NumberOptimizationParameters == 0 || obj.ParentRegion.NumberOptimizationResponses == 0
                return;
            end

            indTorque = findStructName(obj.ParentRegion.OptimizationResponses, "Torque");
            indPhase0 = findStructName(obj.ParentRegion.OptimizationParameters, "Phase0");
            indCurrent = findStructName(obj.ParentRegion.OptimizationParameters, "IApp");
            if ~isempty(indTorque) && ~isempty(indPhase0)
                z_T = z(:, indTorque);
                Lxx(indPhase0, indPhase0, indTorque) = Lxx(indPhase0, indPhase0, indTorque) + z_T'*obj.d2Current_dPhase02(t)*obj.CircuitInterface.getX();
            end
            if ~isempty(indTorque) && ~isempty(indCurrent) && ~isempty(indPhase0)
                z_T = z(:, indTorque);
                crossEntries = z_T'*obj.d2Current_dIAppdPhase0(t)*obj.CircuitInterface.getX();
                Lxx(indPhase0, indCurrent, indTorque) = Lxx(indPhase0, indCurrent, indTorque) + crossEntries;
                Lxx(indCurrent, indPhase0, indTorque) = Lxx(indCurrent, indPhase0, indTorque) + crossEntries;
            end
        end

        function Lxy = function_Lxy(obj, t, x, y, z)
            Lxy = [];
            % Lxy = zeros(obj.ParentRegion.NumberOptimizationParameters, obj.ParentRegion.NumberDOF, obj.ParentRegion.NumberOptimizationResponses);
        end

        function Lyy = function_Lyy(t, x, y, z)
            Lyy = [];
        end
    end
end