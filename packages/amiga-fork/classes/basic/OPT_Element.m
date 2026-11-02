classdef OPT_Element < handle
    properties (Access = public)
        OptimizationResponses, NumberOptimizationResponses;
        OptimizationParameters, NumberOptimizationParameters;
        SaveDirectory;
    end

    methods (Abstract)
        r = responses(obj, t, x, y)
        drdy = partial_r_partial_y(obj, t, x, y)
        drdx = partial_r_partial_x(obj, t, x, y)
        dedx = partial_e_partial_x(obj, t, x, y)
        Lx = function_Lx(obj, t, x, y, z)
        optimizationUpdate(obj, x);
    end

    methods
        function obj = OPT_Element()
            obj.resetOptimizationResponses();
            obj.resetOptimizationParameters();
        end

         function addOptimizationParameter(obj, name, initVal, minVal, maxVal)
            pos = numel(obj.OptimizationParameters) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationParameters.Name)), "This parameter has already been added to the element!");
            end
            assert(minVal~=maxVal, "The minimum and maximum values must be different for parameters that may change!");
            obj.OptimizationParameters(pos).Name = name;
            obj.OptimizationParameters(pos).Init = initVal;
            obj.OptimizationParameters(pos).Min = minVal;
            obj.OptimizationParameters(pos).Max = maxVal;
            obj.NumberOptimizationParameters = obj.NumberOptimizationParameters + 1;
        end

        function addOptimizationResponse(obj, name)
            pos = numel(obj.OptimizationResponses) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationResponses.Name)), "This response has already been added to the element!");
            end
            obj.OptimizationResponses(pos).Name = name;
            obj.NumberOptimizationResponses = obj.NumberOptimizationResponses + 1;
        end

        function resetOptimizationParameters(obj)
            obj.OptimizationParameters = struct('Name', {}, 'Init', {}, 'Min', {}, 'Max', {});
            obj.NumberOptimizationParameters = 0;
        end

        function resetOptimizationResponses(obj)
            obj.OptimizationResponses = struct('Name', {});
            obj.NumberOptimizationResponses = 0;
        end
    end
end