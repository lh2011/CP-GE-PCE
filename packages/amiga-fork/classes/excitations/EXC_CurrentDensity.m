classdef EXC_CurrentDensity < EXC
    properties
        
    end

    methods
        function obj = EXC_CurrentDensity(material, regionNumbers, currentDensityFunction)
            obj@EXC(material, regionNumbers);
            material.addCurrentPostprocessing(obj, obj.RegionNumbers);
            obj.setValue(currentDensityFunction);
        end

        function J = getCurrentDensity(obj, t, x, y, ~, ~)
            if obj.IsSpaceDependent && obj.IsTimeDependent
                spaceFunc = obj.ExcitationFunction(t);
                J = spaceFunc(x,y);
            elseif obj.IsSpaceDependent
                J = obj.ExcitationFunction(x,y);
            elseif obj.IsTimeDependent
                J = ones(size(x))*obj.ExcitationFunction(t);
            else
                J = ones(size(x))*obj.ExcitationFunction;
            end
        end

        function Pv = getLossDensity(obj, t, x, ~, ~, reg)
            Pv = zeros(size(x));
        end
    end
end