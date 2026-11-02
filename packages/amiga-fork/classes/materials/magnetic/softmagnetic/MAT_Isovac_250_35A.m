classdef MAT_Isovac_250_35A < MAT_Magnetic

    methods (Access = public)
        function obj = MAT_Isovac_250_35A(region, number, args)
            arguments
                region 
                number 
                args.LaminationThickness = 0;
            end
            obj@MAT_Magnetic(region, number, "LaminationThickness", args.LaminationThickness);
            obj.MuRel = 6000;
            obj.Sigma = (64.5*1e-8)^-1;
            obj.PlotColor = [0.75, 0.75, 0.75];
            obj.Rho = 7.60;
            obj.IsLinearMAG = false;
            % Values for 50Hz
            obj.Hvalues = [0, 20, 23, 26, 29, 32, 35, 38, 41, 43, 46, 49, 52, 55, 59, 63, 68, 81, 89, 100, 116, 139, 177, 260, 435, 795, 1417, 2263, 3294, 4568, 6131, 7907, 9000, 10000];
            obj.Bvalues = [0, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1731, 1754].*1e-3;
            obj.fitHBspline()
        end
    end
end