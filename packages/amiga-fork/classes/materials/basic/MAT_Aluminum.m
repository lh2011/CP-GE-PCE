classdef MAT_Aluminum < MAT_Magnetic

    methods (Access = public)
        function obj = MAT_Aluminum(region, number)
            obj@MAT_Magnetic(region, number);
            obj.MuRel = 1;
            obj.Sigma = 35e6;
            obj.Rho = 2700;
            obj.PlotColor = "#BCC6CC";
            obj.CostPerKG = 2;
        end
    end
end