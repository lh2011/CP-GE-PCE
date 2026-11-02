classdef MAT_Copper_2nd < MAT_Magnetic_OPT_2nd

    methods (Access = public)
        function obj = MAT_Copper_2nd(region, number)
            obj@MAT_Magnetic_OPT_2nd(region, number);
            obj.MuRel = 1;
            obj.Sigma = 5.96e7;
            obj.Rho = 8960;
            obj.PlotColor = "#E9503E";
            obj.CostPerKG = 10;
        end
    end
end