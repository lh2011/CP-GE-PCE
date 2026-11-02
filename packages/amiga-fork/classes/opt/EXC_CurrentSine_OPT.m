classdef EXC_CurrentSine_OPT < EXC_CurrentSine & OPT_Element
    properties
        dCurrent_dIApp;
        dCurrent_dPhase0;
    end

    methods
        function obj = EXC_CurrentSine_OPT(material, phasesStruct, iPeak, omega, uvwshift, phase0, args)
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
           obj@EXC_CurrentSine(material, phasesStruct, iPeak, omega, uvwshift, phase0, "NumberWindings", args.NumberWindings, "PlotColor", args.PlotColor);
           obj.dCurrent_dIApp = @(t) sin(obj.Omega*t + obj.UVWshift + obj.Phase0);
           obj.dCurrent_dPhase0 = @(t) obj.IApp*cos(obj.Omega*t + obj.UVWshift + obj.Phase0);
        end

        function r = responses(obj, x)
            r = zeros(obj.ParentRegion.NumberOptimizationResponses, 1);
        end

        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberDOF);
        end

        function drdx = partial_r_partial_x(obj, t, x, y)
            drdx = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberOptimizationParameters);
        end

        function Lx = function_Lx(obj, t, x, y, z)
            Lx = zeros(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberOptimizationParameters);
            
            indCurrent = findStructName(obj.ParentRegion.OptimizationParameters, "IApp");
            indTorque = findStructName(obj.ParentRegion.OptimizationResponses, "Torque");
            if ~isempty(indCurrent) && ~isempty(indTorque)
                zT = z(:,indTorque);
                Lx(indTorque, indCurrent) = Lx(indTorque, indCurrent) + obj.dCurrent_dIApp(t)*zT'*obj.CircuitInterface.getX();
            end
            indPhase0 = findStructName(obj.ParentRegion.OptimizationParameters, "Phase0");
            if ~isempty(indPhase0) && ~isempty(indTorque)
                zT = z(:,indTorque);
                Lx(indTorque, indPhase0) = Lx(indTorque, indPhase0) + obj.dCurrent_dPhase0(t)*zT'*obj.CircuitInterface.getX();
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.ParentMaterial.Space.ndof, obj.ParentMaterial.NumberOptimizationParameters);

            indCurrent = findStructName(obj.ParentRegion.OptimizationParameters, "IApp");
            if indCurrent > 0
                dedx(:, indCurrent) = dedx(:, indCurrent) + obj.dCurrent_dIApp(t)*obj.CircuitInterface.getX();
            end
            
            indPhase0 = findStructName(obj.ParentRegion.OptimizationParameters, "Phase0");
            if indPhase0 > 0
                dedx(:, indPhase0) = dedx(:, indPhase0) + obj.dCurrent_dPhase0(t)*obj.CircuitInterface.getX();
            end
        end

        function optimizationUpdate(obj, x)
            if obj.ParentRegion.NumberOptimizationParameters == 0 
                return 
            end

            % TBD CHECK IF THERE IS GEOMETRY-FUNCTION of PARENT-ELEMENT
            % If so: update, adjust also derivatives for change in coil size
            obj.updateExcitation();

            % Update current
            indCurrent = findStructName(obj.ParentRegion.OptimizationParameters, "IApp");
            if indCurrent > 0
                obj.IApp = x(indCurrent);
            end
            % Update phase
            indPhase0 = findStructName(obj.ParentRegion.OptimizationParameters, "Phase0");
            if indPhase0 > 0
                obj.Phase0 = x(indPhase0);
            end
        end
    end
end