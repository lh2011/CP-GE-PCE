classdef MAT_Laminated < MAT_Magnetic
    properties (Access = public) % protected
        f1, f2; % Phase ratio for Iron and Air phases
        BHspline1; % BH characteristic for first phase (pure Iron)
        B1;  % Adjusted B-values for first phase BH spline (pure Iron)
        LamThickness;
    end

    methods (Access = public)

        function fitHBspline(obj)
            % Reverse fitting for nonlinear Iron
            obj.BHspline1 = pchip(obj.H, obj.B1);
            % Adjusted B-values considering air-like lamination
            obj.B = obj.f1*ppval(obj.BHspline1, obj.H) + obj.f2*obj.Mu0*obj.H;

            obj.HBspline = pchip(obj.B, obj.H);
            obj.HBsplineDer = fnder(obj.HBspline);
        end
    end
end