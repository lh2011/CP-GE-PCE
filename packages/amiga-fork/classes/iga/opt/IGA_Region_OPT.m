
classdef IGA_Region_OPT < OPT_System & IGA_Region

    properties (Access=public)
        OptimizationControlPoints;
        GeometryFunction;
        StandardParameters;
        GeometryConstraintFunction; GeometryConstraintName; 
        dCdP;
        % lastx;
    end

    methods
        function obj = IGA_Region_OPT()
            
        end

        function setGeometryFunction(obj, geometryFunction, standardParameters)
            arguments
                obj 
                geometryFunction 
                standardParameters struct = struct("preview", false);
            end
            standardParameters = namedargs2cell(standardParameters);
            obj.GeometryFunction = geometryFunction;
            [~,~,~, obj.StandardParameters] = geometryFunction(standardParameters{:});
        end

        function optimizationUpdate(obj, x)
            % check if model is already up to date
            if ~isempty(obj.lastx) && all(x == obj.lastx)
                return
            end

            if ~isempty(obj.GeometryFunction)
                params = obj.StandardParameters;
                % find which parameters are actually used by the geometry function
                geoNames = fieldnames(obj.StandardParameters);
                for iParam = 1:numel(obj.OptimizationParameters)
                    if any(obj.OptimizationParameters(iParam).Name == geoNames)
                        params.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                    end
                end
                unpackedParams = namedargs2cell(params);
                warning('off', 'geopdes:nrbmultipatch')
                [srf, regions, boundaries] = obj.GeometryFunction(unpackedParams{:});
    
                % Update surface for optimization control points
                ControlPointsUpdate = zeros(size(obj.ControlPoints));
                for iPatch  = 1:obj.SpacesGeo.npatch
                    ind_loc = obj.SpacesGeo.gnum{iPatch};
                    ControlPointsUpdate(ind_loc, 1) = reshape(srf(iPatch).coefs(1, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsUpdate(ind_loc, 2) = reshape(srf(iPatch).coefs(2, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsUpdate(ind_loc, 3) = reshape(srf(iPatch).coefs(3, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsUpdate(ind_loc, 4) = reshape(srf(iPatch).coefs(4, :, :, :), [], 1);
                end
                % Apply control point offsets
                for iCtrlPntOffsets = 1:numel(obj.OptimizationControlPoints)
                    paramNr = findStructName(obj.OptimizationParameters, obj.OptimizationControlPoints(iCtrlPntOffsets).Name);
                    dim = obj.OptimizationControlPoints(iCtrlPntOffsets).Dimension;
                    for iCtrlPnt = 1:numel(obj.OptimizationControlPoints(iCtrlPntOffsets).Numbers)
                        CtrlPntNr = obj.OptimizationControlPoints(iCtrlPntOffsets).Numbers(iCtrlPnt);
                        if dim == "x"
                            ControlPointsUpdate(CtrlPntNr, 1) = ControlPointsUpdate(CtrlPntNr, 1) + x(paramNr);
                        elseif dim == "y"
                            ControlPointsUpdate(CtrlPntNr, 2) = ControlPointsUpdate(CtrlPntNr, 2) + x(paramNr);
                        elseif dim == "r"
                            angle = obj.OptimizationControlPoints(iCtrlPntOffsets).Angles(iCtrlPnt);
                            ControlPointsUpdate(CtrlPntNr, 1) = ControlPointsUpdate(CtrlPntNr, 1) + cos(angle)*x(paramNr);
                            ControlPointsUpdate(CtrlPntNr, 2) = ControlPointsUpdate(CtrlPntNr, 2) + sin(angle)*x(paramNr);
                        end
                    end
                end
                % check if evaluated first time
                if isempty(obj.dCdP)
                    obj.updatedCdP(x);
                end
            
                % only update meshes and spaces if there are geometry changes
                if ~all(ControlPointsUpdate == obj.ControlPoints, 'all')
                    % update surface with offsets
                    for iPatch  = 1:obj.SpacesGeo.npatch
                        ind_loc = obj.SpacesGeo.gnum{iPatch};
                        srf(iPatch).coefs(1, :, :, :) = reshape(ControlPointsUpdate(ind_loc, 1).*ControlPointsUpdate(ind_loc, 4), srf(iPatch).number);
                        srf(iPatch).coefs(2, :, :, :) = reshape(ControlPointsUpdate(ind_loc, 2).*ControlPointsUpdate(ind_loc, 4), srf(iPatch).number);
                        srf(iPatch).coefs(3, :, :, :) = reshape(ControlPointsUpdate(ind_loc, 3).*ControlPointsUpdate(ind_loc, 4), srf(iPatch).number);
                        srf(iPatch).coefs(4, :, :, :) = reshape(ControlPointsUpdate(ind_loc, 4), srf(iPatch).number);
                    end
                    
                    obj.importSurface(srf, regions, boundaries);
                    obj.updatedCdP(x);
                end
            end

            for iRegion = 1:numel(obj.Regions)
                obj.Regions(iRegion).Material.optimizationUpdate(x);
            end
            obj.lastx = x;
        end
        
        function r = responses(obj, t, x, y)
            obj.optimizationUpdate(x)
            % get responses from elements (materials)
            r = responses@OPT_System(obj, t, x, y);

            indSmoothness = findStructName(obj.OptimizationResponses, "SurfaceSmoothness");
            if ~isempty(indSmoothness)
                r(indSmoothness) = obj.surfaceSmoothness(x);
            end

            if ~isempty(obj.GeometryConstraintName)
                xStruct = obj.xToStruct(x);
                unpackedStruct = namedargs2cell(xStruct);
                c = obj.GeometryConstraintFunction(unpackedStruct{:});
                indGeometryConstraints = contains(vertcat(obj.OptimizationResponses.Name), obj.GeometryConstraintName);
                r(indGeometryConstraints) = c;
            end
        end
        
        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = sparse(obj.NumberOptimizationResponses, obj.NumberDOF);
        end

        function drdx = partial_r_partial_x(obj, t, x, y)
            obj.optimizationUpdate(x)

            % get responses from elements (materials)
            drdx = partial_r_partial_x@OPT_System(obj, t, x, y);

            indSmoothness = findStructName(obj.OptimizationResponses, "SurfaceSmoothness");
            if ~isempty(indSmoothness)
                [~, dSdx] = obj.surfaceSmoothness(x);
                drdx(indSmoothness, :) = dSdx;
            end
            % Geometry Constraints: numerical evaluation of gradient
            if ~isempty(obj.GeometryConstraintName)
                indGeometryConstraints = contains(vertcat(obj.OptimizationResponses.Name), obj.GeometryConstraintName);
                dcdx = obj.geometryConstraintGradient(x);
                drdx(indGeometryConstraints, :) = dcdx;
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            obj.optimizationUpdate(x);
            % get responses from elements (materials)
            dedx = partial_e_partial_x@OPT_System(obj, t, x, y);
        end

        function Lx = function_Lx(obj, t, x, y, z)
            obj.optimizationUpdate(x);
            % get responses from elements (materials)
            Lx = function_Lx@OPT_System(obj, t, x, y, z);
        end
        
        function [p, dSdx, d2Sdx2] = surfaceSmoothness(obj, x)
            r = zeros(numel(obj.OptimizationControlPoints), 1);
            angles = zeros(numel(obj.OptimizationControlPoints), 1);
            r_ind = zeros(numel(obj.OptimizationControlPoints), 1);
            
            % find all angles and offsets first
            for iCtrlPnt = 1:numel(obj.OptimizationControlPoints)
                CtrlPntNr = obj.OptimizationControlPoints(iCtrlPnt).Numbers(1);
                angles(iCtrlPnt) = atan2(obj.ControlPoints(CtrlPntNr, 2), obj.ControlPoints(CtrlPntNr, 1));
                paramNr = find(obj.OptimizationControlPoints(iCtrlPnt).Name == vertcat(obj.OptimizationParameters.Name));
                r_ind(iCtrlPnt) = paramNr;
                r(iCtrlPnt) = x(paramNr);
            end
            [~, sort_ind] = sort(angles);
            deltar = diff(angles(sort_ind));
            % Penalty matrix which containts angle differences
            K = zeros(numel(r), numel(r));
            for iCtrlPnt = 1:numel(deltar)
                K(iCtrlPnt, iCtrlPnt) = K(iCtrlPnt, iCtrlPnt) + 1./deltar(iCtrlPnt);
                K(iCtrlPnt+1, iCtrlPnt) = K(iCtrlPnt+1, iCtrlPnt) - 1./deltar(iCtrlPnt);
                K(iCtrlPnt, iCtrlPnt+1) = K(iCtrlPnt, iCtrlPnt+1) - 1./deltar(iCtrlPnt);
                K(iCtrlPnt+1, iCtrlPnt+1) = K(iCtrlPnt+1, iCtrlPnt+1) + 1./deltar(iCtrlPnt);
            end
            p = r(sort_ind)'*K*r(sort_ind);

            dSdx = zeros(obj.NumberOptimizationParameters, 1);
            dSdx(r_ind(sort_ind)) = 2*K*r(sort_ind);
            d2Sdx2(r_ind(sort_ind), r_ind(sort_ind)) = 2*K;
        end

        function dcdx = geometryConstraintGradient(obj, x)
            xStruct = obj.xToStruct(x);
            unpackedStruct = namedargs2cell(xStruct);
            c0 = obj.GeometryConstraintFunction(unpackedStruct{:});
            
            dcdx = zeros(numel(c0), obj.NumberOptimizationParameters);
            for iParam = 1:obj.NumberOptimizationParameters
                x1 = x;
                Step = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-6;
                x1(iParam) = x1(iParam) + Step;
                xStruct1 = obj.xToStruct(x1);
                unpackedStruct1 = namedargs2cell(xStruct1);
                c1 = obj.GeometryConstraintFunction(unpackedStruct1{:});
                dcdx(:, iParam) = (c1-c0)/Step;
            end
        end

        function updatedCdP(obj, x)
            % update only if there is a geometry function set (aka it shall be optimized)
            if isempty(obj.GeometryFunction)
                obj.dCdP = [];
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
                    steps(usedParams) = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-6;
                    indUsedParams = [indUsedParams, iParam];
                end
            end

            obj.dCdP = zeros([size(obj.ControlPoints), obj.NumberOptimizationParameters]);
            if usedParams > 0
                obj.dCdP(:,:,indUsedParams) = getdCdP(obj.GeometryFunction, obj.StandardParameters, checkParams, obj.SpacesGeo, steps);
            end
            % check optimization control points
            for iCtrlPntOffsets = 1:numel(obj.OptimizationControlPoints)
                paramNr = findStructName(obj.OptimizationParameters, obj.OptimizationControlPoints(iCtrlPntOffsets).Name);
                dim = obj.OptimizationControlPoints(iCtrlPntOffsets).Dimension;
                for iCtrlPnt = 1:numel(obj.OptimizationControlPoints(iCtrlPntOffsets).Numbers)
                    CtrlPntNr = obj.OptimizationControlPoints(iCtrlPntOffsets).Numbers(iCtrlPnt);
                    if dim == "x"
                        obj.dCdP(CtrlPntNr, 1, paramNr) = 1;
                    elseif dim == "y"
                        obj.dCdP(CtrlPntNr, 2, paramNr) = 1;
                    elseif dim == "r"
                        angle = obj.OptimizationControlPoints(iCtrlPntOffsets).Angles(iCtrlPnt);
                        obj.dCdP(CtrlPntNr, 1, paramNr) = cos(angle);
                        obj.dCdP(CtrlPntNr, 2, paramNr) = sin(angle);
                    end
                end
            end
        end

        function addOptimizationResponse(obj, name)
            addOptimizationResponse@OPT_Element(obj, name);
            % pass response name to the materials
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.addOptimizationResponse(name);
                end
            end
            obj.lastx = [];
            obj.dCdP = [];
        end

        function resetOptimizationResponses(obj)
            resetOptimizationResponses@OPT_Element(obj);
            % pass to the materials
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.resetOptimizationResponses();
                end
            end
        end

        function resetOptimizationParameters(obj)
            resetOptimizationParameters@OPT_Element(obj);
            % pass to the materials
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.resetOptimizationParameters();
                end
            end
            obj.OptimizationControlPoints = struct('Name', {}, 'Numbers', {}, 'Dimension', {});
        end

        function addGeometryConstraintFunction(obj, name, func)
            obj.GeometryConstraintFunction = func;
            ndim = numel(func());
            obj.GeometryConstraintName = name;
            for iDim = 1:ndim
                obj.addOptimizationResponse(name+num2str(iDim));
            end
        end

        % Functions for adding responses to the element
        function addOptimizationParameter(obj, name, initVal, minVal, maxVal)
            addOptimizationParameter@OPT_Element(obj, name, initVal, minVal, maxVal);
            % pass parameter name to the materials
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.addOptimizationParameter(name, initVal, minVal, maxVal);
                end
            end
            obj.lastx = [];
            obj.dCdP = [];
        end

        % dim: "x", "y", "r"
        function addOptimizationControlPointOffset(obj, name, numbers, dim, initVal, minVal, maxVal)
            assert(any(dim==["x", "y", "r"]), "The dimension mus be either x, y or r!");
            obj.addOptimizationParameter(name, initVal, minVal, maxVal);
            pos = numel(obj.OptimizationControlPoints) + 1;
            obj.OptimizationControlPoints(pos).Name = name;
            obj.OptimizationControlPoints(pos).Numbers = numbers;
            obj.OptimizationControlPoints(pos).Dimension = dim;
            % Angle for radial offset fixed from beginning
            obj.OptimizationControlPoints(pos).Angles = atan2(obj.ControlPoints(numbers, 2), obj.ControlPoints(numbers, 1));
        end

        function plotOptimizationControlPoints(obj)
            numbers = [];
            for iCtrlPts = 1:numel(obj.OptimizationControlPoints)
                numbers = [numbers, reshape(obj.OptimizationControlPoints(iCtrlPts).Numbers, 1, [])];
            end
            scatter(obj.ControlPoints(numbers, 1), obj.ControlPoints(numbers, 2), "filled", "red")
            text(obj.ControlPoints(numbers, 1), obj.ControlPoints(numbers, 2), string(numbers), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
        end

        function plotFunction(obj, x, optimValues, state)
            obj.optimizationUpdate(x);
            % set(gcf, "Name", "Geometry")
            obj.plotGeometry("FigureName", "Optimization Plot Function", "FocusFigure", false);
            nrbexport(obj.Surface, obj.SaveDirectory + "/" + obj.Name+"Opt"+string(optimValues.iteration)+".txt")
        end

        function xStruct = xToStruct(obj, x)
            xStruct = obj.StandardParameters;
            geoNames = fieldnames(obj.StandardParameters);
            for iParam = 1:obj.NumberOptimizationParameters
                if any(obj.OptimizationParameters(iParam).Name == geoNames)
                    xStruct.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                end
            end
        end
    end
end

 