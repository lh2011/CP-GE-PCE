classdef MAT_IronAnalytic2 < MAT_Magnetic
    properties
        Hs = 100;
        Bs = 2;
    end

    methods (Access = public)
        function obj = MAT_IronAnalytic2()
            obj.Mur = 2000;
            obj.IsLinearMAG = false;
            obj.Sigma = 1.03e7;
        end

        function mu = getMuNonlinear(obj, B)
            mu = 1./obj.getNuNonlinear(B);
        end

        function nu = getNuNonlinear(obj, B)
            nu = obj.Hs/obj.Bs^2.*B.^2 + 1e-3;
            % nu = ones(size(B))*1000;
        end

        function nuPrime = getNuPrimeNonlinear(obj, B)
            nuPrime = 2*B.*obj.Hs/obj.Bs^2.*ones(size(B));
            % nuPrime = ones(size(B))*0;
        end

        function d2nudB2 = getd2nudB2(obj, B)
            d2nudB2 = 2*obj.Hs/obj.Bs^2.*ones(size(B));
            % d2nudB2 = zeros(size(B));
        end
    end
end