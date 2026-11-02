classdef OPT_Element_2nd < OPT_Element
    properties (Access = public)
        
    end

    methods (Abstract)
        % Lxx = zeros    [#nOptParam, #nOptParam, #nOptResponses] 
        Lxx = function_Lxx(obj, t, x, y, z)
        % Lxy = sptensor [#nOptParam, #nDOF, #nOptResponses] 
        Lxy = function_Lxy(obj, t, x, y, z)
        % Lyy = sptensor [#nDOF, #nDOF, #nOptResponses] 
        Lyy = function_Lyy(obj, t, x, y, z)
    end
end