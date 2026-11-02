classdef MAT_BrauerIron < MAT_Magnetic
    properties
        k1 = 0.3774
        k2 = 2.970
        k3 = 388.33
        Bclip = 10;
        muClip;
%         Hclip;
    end

    methods (Access = public)
        function obj = MAT_BrauerIron()
            obj.Mur = 1200;
            obj.IsLinearMAG = false;
            obj.Sigma = 1.03e7;
            obj.Rho = 7.874;
            obj.correctBclip();
            
        end
        % Calculate the magetic flux density, where a further decreasing of
        % the differential permeability is not physical (Bclip)
        % permeabilities higher than Bclip will increase linearly with mu0
        function correctBclip(obj)
            b = 0:0.001:3;
            h = b.*obj.getNuNonlinear(b);
            mudiff = gradient(b, h);
            obj.Bclip = b(find(mudiff<obj.Mu0, 1));
            obj.muClip = obj.getMuNonlinear(obj.Bclip);
            obj.B = linspace(0,obj.Bclip, 50);
            obj.H = obj.B.*obj.getNuNonlinear(obj.B);
        end

        function nu = getNuNonlinear(obj, B)
            nu = obj.k1*exp(obj.k2.*B.^2) + obj.k3;
            % correction for too high unphysical values (mu continues linearly)
            nu(B>obj.Bclip) = (B(B>obj.Bclip)-obj.Bclip.*(1-obj.Mu0./obj.muClip))./(B(B>obj.Bclip).*obj.Mu0);
        end
        function mu = getMuNonlinear(obj, B)
            mu = 1./obj.getNuNonlinear(B);
        end

        function d2nudB2 = getd2nudB2(obj, B)
            d2nudB2 = 2*obj.k1.*obj.k2.*exp(obj.k2.*B.^2) + 4*B.^2.*obj.k1.*obj.k2.^2.*exp(obj.k2.*B.^2);

            B2correct = B(B>obj.B(end));
            d2nudB2(B>obj.B(end)) = -2*(obj.B(end)-obj.Mu0*obj.H(end))/obj.Mu0./B2correct.^3;         
        end

        function nuPrime = getNuPrimeNonlinear(obj, B)
            nuPrime = 2*B.*obj.k1.*obj.k2.*exp(obj.k2.*B.^2);

            B2correct = B(B>obj.B(end));
            nuPrime(B>obj.B(end)) = (obj.B(end)-obj.Mu0*obj.H(end))/obj.Mu0./B2correct.^2; 
        end
    end
end