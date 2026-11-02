classdef IGA_Region_OPT_2nd < IGA_Region_OPT & OPT_System_2nd

    properties        
        d2CdP2;
        lastx_2nd;
    end

    methods
        function optimizationUpdate2nd(obj, x)
            % Check if update is necessary!
            if isempty(obj.lastx_2nd) || any(x~=obj.lastx_2nd)
                obj.updated2CdP2(x);
                obj.lastx_2nd = x;
            end
        end
        
        function Lxx = function_Lxx(obj, t, x, y, z)
            obj.optimizationUpdate2nd(x);
            Lxx = function_Lxx@OPT_System_2nd(obj, t, x, y, z);
             
            indSmoothness = findStructName(obj.OptimizationResponses, "SurfaceSmoothness");
            if ~isempty(indSmoothness)
                [~,~,d2Sdx2] = obj.surfaceSmoothness(x);
                Lxx(:,:,indSmoothness) = d2Sdx2;
            end

            if ~isempty(obj.GeometryConstraintName)
                indGeometryConstraints = find(contains(vertcat(obj.OptimizationResponses.Name), obj.GeometryConstraintName));
                d2cdx2 = obj.geometryConstraintGradient_2nd(x);
                Lxx(:, :, indGeometryConstraints') = d2cdx2;
            end
        end

        function d2cdx2 = geometryConstraintGradient_2nd(obj, x)
            dcdx = obj.geometryConstraintGradient(x);
            d2cdx2 = zeros(obj.NumberOptimizationParameters, obj.NumberOptimizationParameters, size(dcdx, 1));
            for iParam = 1:obj.NumberOptimizationParameters
                x1 = x;
                Step = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-6;
                x1(iParam) = x1(iParam) + Step;
                dcdx1 = obj.geometryConstraintGradient(x1);
                d2cdx2(iParam,:,:) = (dcdx1-dcdx)'/Step;
            end
        end

        function updated2CdP2(obj, x)
            % update only if there is a geometry function set (aka it shall be optimized)
            if isempty(obj.GeometryFunction)
                obj.d2CdP2 = [];
                return
            end

            % find which parameters are actually used by the geometry function
            geoNames = fieldnames(obj.StandardParameters);
            usedParams = 0;
            indUsedParams = [];
            for iParam = 1:numel(obj.OptimizationParameters)
                if any(obj.OptimizationParameters(iParam).Name == geoNames)
                    usedParams = usedParams+1;
                    checkParams.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                    % make numerical step size relative to parameter ranges
                    steps(usedParams) = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-5;
                    indUsedParams = [indUsedParams, iParam];
                end
            end

            obj.d2CdP2 = zeros([size(obj.ControlPoints), obj.NumberOptimizationParameters, obj.NumberOptimizationParameters]);
            mytimer("d2CdP2", "start");
            obj.d2CdP2(:,:,indUsedParams,indUsedParams) = getd2CdP2(obj.GeometryFunction, obj.StandardParameters, checkParams, obj.SpacesGeo, steps);
            mytimer("d2CdP2", "stop");
        end
    end
end

 