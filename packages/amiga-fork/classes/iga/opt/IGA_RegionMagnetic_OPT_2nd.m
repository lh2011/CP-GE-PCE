classdef IGA_RegionMagnetic_OPT_2nd <  IGA_RegionMagnetic_OPT & IGA_Region_OPT_2nd

    properties
        
    end

    methods
        function obj = IGA_RegionMagnetic_OPT_2nd(subdivisions, degrees, length, args)
            arguments
                subdivisions (1, 2) double  = [5, 5];
                degrees (1, 2) double = [2, 2];
                length (1, 1) double = 1;
                args.degreeQuadrature = degrees + [1, 1];
            end
            obj@IGA_RegionMagnetic_OPT(subdivisions, degrees, length, "degreeQuadrature",args.degreeQuadrature);
        end
    end
end

