classdef IGA_RegionMagnetic_OPT <  IGA_RegionMagnetic & IGA_Region_OPT

    properties
        OptimizationEvaluationPoints = struct("Name", {}, "Numbers", {});
    end

    methods
        function obj = IGA_RegionMagnetic_OPT(subdivisions, degrees, length, args)
            arguments
                subdivisions (1, 2) double  = [5, 5];
                degrees (1, 2) double = [2, 2];
                length (1, 1) double = 1;
                args.degreeQuadrature = degrees + [1, 1];
            end
            obj@IGA_RegionMagnetic(subdivisions, degrees, length, "degreeQuadrature", args.degreeQuadrature);
        end

        function r = responses(obj, t, x, y)
            r = responses@IGA_Region_OPT(obj, t, x, y);

            % B magnitude of evaluation points
            for iEvalPoint = 1:numel(obj.OptimizationEvaluationPoints)
                name = obj.OptimizationEvaluationPoints(iEvalPoint).Name;
                numbers = obj.OptimizationEvaluationPoints(iEvalPoint).Numbers;
                ind_beval = find (vertcat(obj.OptimizationResponses.Name) == name);
                if ~isempty(ind_beval)
                    r(ind_beval) = obj.BMagEvalPoints(y, numbers);
                end
            end
        end

        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = partial_r_partial_y@IGA_Region_OPT(obj, t, x, y);

            for iEvalPoint = 1:numel(obj.OptimizationEvaluationPoints)
                name = obj.OptimizationEvaluationPoints(iEvalPoint).Name;
                numbers = obj.OptimizationEvaluationPoints(iEvalPoint).Numbers;
                ind_beval = find (vertcat(obj.OptimizationResponses.Name) == name);
                if ~isempty(ind_beval)
                    [~, Bmagder] = obj.BMagEvalPoints(y, numbers);
                    drdy(ind_beval, :) = Bmagder;
                end
            end
        end

        function addOptimizationResponseBeval(obj, name, numbers)
            if ~(numel(intersect(numbers, 1:size(obj.EvaluationPoints, 2))) == numel(numbers))
                error("You need to define the optimization points as evaluation points first");
            end
            obj.addOptimizationResponse(name);
            pos = numel(obj.OptimizationEvaluationPoints) + 1;
            obj.OptimizationEvaluationPoints(pos).Name = name;
            obj.OptimizationEvaluationPoints(pos).Numbers = numbers;
        end

        function [Bmag, Bmagder] = BMagEvalPoints(obj, y, numbers)
            GapShapeGradX = reshape(sparse(double(obj.EvaluationShapeValues{2}(numbers,1,:))), numel(numbers), []);
            GapShapeGradY = reshape(sparse(double(obj.EvaluationShapeValues{2}(numbers,2,:))), numel(numbers), []);

            BmagValues = ((GapShapeGradX*y).^2 + (GapShapeGradY*y).^2).^0.5;
            Bmag = mean(BmagValues);
            if nargout == 2
                Ni_x_ui = GapShapeGradX*y;
                Ni_y_ui = GapShapeGradY*y;
                Bmagder = sum((GapShapeGradX.*Ni_x_ui +  GapShapeGradY.*Ni_y_ui)./BmagValues, 1)./numel(numbers);
            end
        end
    end
end

