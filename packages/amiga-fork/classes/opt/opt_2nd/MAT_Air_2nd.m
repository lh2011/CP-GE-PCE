classdef MAT_Air_2nd < MAT_Magnetic_OPT_2nd

    methods (Access = public)
        function obj = MAT_Air_2nd(region, number)
            obj@MAT_Magnetic_OPT_2nd(region, number);
            %https://theengineeringmindset.com/properties-of-air-at-atmospheric-pressure/
            obj.MuRel = 1;
            obj.Sigma = 0;
            obj.Rho = 1.2047;
            obj.PlotColor = "#009CDA";

            obj.CostPerKG = 0;
        end
    end
end