classdef MAT_Magnet < MAT_Magnetic
    properties (Access = public)
        Br;
        Angle;
    end
    methods (Access = public)
        function obj = MAT_Magnet(region, number, remanenceDirection, args)
            arguments
                region;
                number;
                remanenceDirection;
                args.Remanence = 1;
                args.RelativePermeability = 1.05;
            end
            obj@MAT_Magnetic(region, number);
            obj.MuRel = args.RelativePermeability;
            obj.Br = args.Remanence;
            if isnumeric(remanenceDirection)
                obj.Angle = @(x,y) ones(size(x))*remanenceDirection;
            elseif isa(remanenceDirection, 'function_handle') && nargin(remanenceDirection) == 2
                obj.Angle = remanenceDirection;
            else
                error("MAT_Magnet:The remanence direction must be either a scalar or a function depending on x and y");
            end
            % Standard values
            obj.CostPerKG = 50;
            obj.MuRel = 1.05;
            obj.Sigma = 1/(1.4e-6);
            obj.Rho = 7600;
            obj.PlotColor = "#00715E";
        end
        
        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.Space.ndof, 1);
            if nargout == 2
                Flin = obj.Mesh.length*int_gradv_Br_bot(obj.Space, obj.Mesh, obj);
            end
        end
        
        function [Bx, By] = getRemanence(obj, x, y)
            Bx = cos(obj.Angle(x,y))*obj.Br;
            By = sin(obj.Angle(x,y))*obj.Br;
        end

        function [Hx, Hy] = getMagneticFieldStrength(obj, Bx, By, x, y)
            alpha = obj.Angle(x,y);
            nu = obj.getReluctivity((Bx.^2+By.^2).^0.5);
            Hx = nu.*(Bx - obj.Br*cos(alpha));
            Hy = nu.*(By - obj.Br*sin(alpha));
        end

        function E = getEnergyDensity(obj, Bx, By, x, y)
            alpha = obj.Angle(x,y);
            % rotation in reference system
            Balpha = Bx.*cos(alpha) + By.*sin(alpha);
            Bbeta = -Bx.*sin(alpha) + By.*cos(alpha);

            % linear for now
            Brem = obj.Br;

            B = (Bx.^2+By.^2).^0.5;
            nu = obj.getReluctivity(B); 
            % int(B.dH) from Br to Balpha. see https://www.femm.info/wiki/PMEnergy#ref3
            E = 0.5.*nu.*(Balpha - Brem).^2 + 0.5.*nu.*Bbeta.^2; 
        end

        function C = getCoenergyDensity(obj, Bx, By, x, y)
            alpha = obj.Angle(x,y);
            % rotation in reference system
            Balpha = Bx.*cos(alpha) + By.*sin(alpha);
            Bbeta = -Bx.*sin(alpha) + By.*cos(alpha);

            % linear for now
            B = (Bx.^2+By.^2).^0.5;
            nu = obj.getReluctivity(B); 
            % int(H.dB) from -Hc to H. see https://www.femm.info/wiki/PMEnergy#ref3
            C = 0.5.*nu.*Balpha.^2 + 0.5.*nu.*Bbeta.^2; 
        end
    end
end