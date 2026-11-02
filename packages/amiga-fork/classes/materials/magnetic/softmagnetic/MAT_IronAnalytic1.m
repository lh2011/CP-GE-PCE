classdef MAT_IronAnalytic1 < MAT_Magnetic
    properties
        J0 = 1.5;
        Bclip = 10;
        muClip;
%         Hclip;
    end

    methods (Access = public)
        function obj = MAT_IronAnalytic1()
            obj.Mur = 2000;
            obj.IsLinearMAG = false;
            obj.Sigma = 1.03e7;
            obj.Rho = 7.874;
            obj.H = linspace(0, 1e5, 1e4);
            obj.B = obj.Mu0*obj.H + (2*obj.J0)./pi* atan(pi*(obj.Mur-1)*obj.Mu0.*obj.H/(2*obj.J0));
            obj.fitHBspline();
        end
    end
end