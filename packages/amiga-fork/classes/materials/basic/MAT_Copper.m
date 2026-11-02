classdef MAT_Copper < MAT_Magnetic_OPT

    methods (Access = public)
        function obj = MAT_Copper(region, number)
            obj@MAT_Magnetic_OPT(region, number);
            obj.MuRel = 1;
            obj.Sigma = 5.96e7;
            obj.Rho = 8960;
            obj.PlotColor = "#E9503E";
            obj.CostPerKG = 10;
        end
    end
end