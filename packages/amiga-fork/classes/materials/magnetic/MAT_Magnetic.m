classdef MAT_Magnetic < MAT & EQS_Element
    properties
        IsLinearMAG = true;
        Mu0 = 4*pi*1e-7;
        Sigma = 0;
        MuRel = 1;
        Bvalues, Hvalues;
        HBspline;
        HBsplineDer;
        HBsplineDer2;
        HBsplineInt;
        IsLaminatedZ;
        LaminationThickness;
        FillFactor;
        CurrentPostprocessing;
    end

    methods (Access = public)
        function obj = MAT_Magnetic(region, regionNumbers, args)
            arguments
                region;
                regionNumbers;
                args.FEMMname = [];
                args.RelativePermeability = 1;
                args.Conductivity = 0;
                args.Linear (1,1) logical = true;
                args.LaminationThickness = 0;
                args.FillFactor = 1;
                args.PlotColor = [0.7,0.7,0.7];
            end
            obj.Region = region;
            obj.RegionNumbers = regionNumbers;
            obj.IsLaminatedZ = args.LaminationThickness > 0;
            obj.FillFactor = args.FillFactor;
            obj.PlotColor = args.PlotColor;
            obj.setLinear(args.Linear);
            obj.setRelativePermeability(args.RelativePermeability);
            obj.setConductivity(args.Conductivity);
            obj.LaminationThickness = args.LaminationThickness;
            % condensing MESH and SPACE coming from the geometry class
            elementNumbers = region.getPatchNumbers(regionNumbers);
            obj.Mesh = MESH_2D(region.MeshesEval(elementNumbers), "length", region.Length);
            obj.Space = SPACE_2D(region.SpacesEval(elementNumbers), region.Spaces.gnum(elementNumbers), region.Spaces.ncomp);
            obj.NumberDOF = 0;                  % material has no intrinsic dofs but writes to the external space object
            region.setMaterial(obj, regionNumbers);    % assign material to region

            if ~isempty(args.FEMMname)
                obj.loadFEMMmaterial(args.FEMMname);
            end
            obj.resetCurrentPostprocessing();
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            if nargout == 2
                if obj.IsLaminatedZ
                    coeff = 1/12*obj.getConductivity()*obj.LaminationThickness^2;
                    Mlin = obj.Mesh.length*int_gradu_gradv(obj.Space, obj.Space, obj.Mesh, coeff);
                else
                    Mlin = obj.Mesh.length*int_u_v(obj.Space, obj.Space, obj.Mesh, obj.getConductivity());
                end
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            % if obj.IsLinearMAG
                Knonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            % else
                % Knonlin = - obj.Mesh.length*int_gradu_nu_gradv(obj.Space, obj.Space, obj.Mesh, y, obj);
            % end

            if nargout == 2
                if obj.IsLinearMAG
                    nu = obj.getReluctivity();
                    Klin = - obj.Mesh.length*int_gradu_gradv(obj.Space, obj.Space, obj.Mesh, nu);
                else
                    Klin = sparse(obj.Space.ndof, obj.Space.ndof);
                end
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            if obj.IsLinearMAG
                Fnonlin = zeros(obj.Space.ndof, 1);
            else
                Fnonlin = - obj.Mesh.length*int_gradu_nu_gradv_times_u(obj.Space, obj.Space, obj.Mesh, y, obj);
            end
            Flin = zeros(obj.Space.ndof, 1);
            for iExc = 1:numel(obj.Excitations)
                Exc = obj.Excitations(iExc).Excitation;
                if nargout == 2
                    [Fnonlin_exc, Flin_exc] = Exc.getEXCvalues(t);
                    Fnonlin(Exc.LocalIndices) = Fnonlin(Exc.LocalIndices) + Fnonlin_exc;
                    Flin(Exc.LocalIndices) = Flin(Exc.LocalIndices) + Flin_exc;
                else
                    Fnonlin_exc = Exc.getEXCvalues(t);
                    Fnonlin(Exc.LocalIndices) = Fnonlin(Exc.LocalIndices) + Fnonlin_exc;
                end
            end
        end

        function J = JacobianMatrix(obj, t, y)
            if obj.IsLinearMAG
                J = sparse(obj.Space.ndof, obj.Space.ndof);
            else
                J = - obj.Mesh.length*int_d_dy_gradu_u_nu_gradv(obj.Space, obj.Space, obj.Mesh, y, obj);
            end
        end

        function resetCurrentPostprocessing(obj)
            obj.CurrentPostprocessing = struct('Element', {}, 'Regions', {});
        end

        function addCurrentPostprocessing(obj, el, regions)
            obj.CurrentPostprocessing(end+1).Element = el;
            obj.CurrentPostprocessing(end).Regions = regions;
        end

        function postprocess(obj, t, y, yprime)
            obj.Solution.Time = t;
            obj.Solution.Values = y;
            obj.Solution.ValuesDer = yprime;
        end

        function loadFEMMmaterial(obj, femmName)
            obj.Name = femmName;
            materialStruct = parseFEMMmagnetic(femmName);
            if isempty(materialStruct)
                error('MAT_Magnetic: FEMM Material not found');
            end
            obj.setConductivity(materialStruct.Sigma*1e6);
            obj.setRelativePermeability(materialStruct.Mu_x);
            if isempty(materialStruct.BH)
                obj.IsLinearMAG = true;
            else
                obj.IsLinearMAG = false;
                obj.Bvalues = materialStruct.BH(:,1)';
                obj.Hvalues = materialStruct.BH(:,2)';
                obj.fitHBspline(true);
            end
        end

        function setLinear(obj, trueFalse)
            arguments
                obj
                trueFalse (1,1) logical
            end
            obj.IsLinearMAG = trueFalse;
        end

        function mu = getPermeability(obj, B)
            arguments
                obj 
                B = 0;
            end
            if obj.IsLinearMAG
                mu = obj.MuRel*obj.Mu0*ones(size(B));
            else
                B(B<=1e-8) = 1e-8; % avoid zero division
                mu = B./ppval(obj.HBspline, B);
                % correct B values higher than measurement (assume behavior like vacuum)
                mu(B>obj.Bvalues(end)) = B(B>obj.Bvalues(end))./((B(B>obj.Bvalues(end))-obj.Bvalues(end))./obj.Mu0 + obj.Hvalues(end));
            end
        end

        function nu = getReluctivity(obj, B)
            arguments
                obj 
                B = 0;
            end
            if obj.IsLinearMAG
                nu = 1/(obj.MuRel*obj.Mu0)*ones(size(B));
            else
                nu = 1./obj.getPermeability(B);
            end
        end

        function sig = getConductivity(obj)
            sig = obj.Sigma;
        end

        % default function for isotropic material
        function [Bx, By] = getRemanence(obj, x, y)
            Bx = zeros(size(x));
            By = zeros(size(y));
        end

        % default function for isotropic material
        function [Hx, Hy] = getMagneticFieldStrength(obj, Bx, By, x, y)
            B = (Bx.^2+Bx.^2).^0.5;
            nu = obj.getReluctivity(B);
            Hx = nu.*Bx;
            Hy = nu.*By;
        end

        % default function for isotropic material
        function E = getEnergyDensity(obj, Bx, By, ~, ~)
            B = (Bx.^2 + By.^2).^0.5;
            if obj.IsLinearMAG
                E = 1./(2*obj.getPermeability(B)).*B.^2;
            else
                E = ppval(obj.HBsplineInt, B);
            end
        end

        % default function for isotropic material
        function C = getCoenergyDensity(obj, Bx, By, ~, ~)
            if obj.IsLinearMAG
                C = 1./(2*obj.getPermeability(B)).*B.^2;
            else
                C = B.^2./obj.getPermeability(B) - obj.getEnergyDensity(Bx, By);
            end
        end

        function D = getDemagnetization(obj, Bx, By, x, y)
            D = zeros(size(x));
        end

        function Jz = getCurrentDensity(obj, t, AzPrime, reg, x, y)
            % default: plain massive solid
            if obj.IsLaminatedZ || obj.getConductivity() == 0
                Jz = zeros(size(x));
            else
                Jz = obj.Sigma.*AzPrime;
            end
            % in case of current density / stranded / solid conductor: overwrite
            for iEl = 1:numel(obj.CurrentPostprocessing)
                el = obj.CurrentPostprocessing(iEl).Element;
                if ismember(reg, obj.CurrentPostprocessing(iEl).Regions)
                    Jz = el.getCurrentDensity(t, x, y, AzPrime, reg);
                end
            end
        end

        function Jz = getLossDensity(obj, t, AzPrime, reg, x, y)
            if obj.IsLaminatedZ || obj.getConductivity() == 0
                warning("TBD calculate correctly!")
                Jz = zeros(size(x));
            else
                % default: massive solid
                Jz = obj.Sigma.*AzPrime.^2;
            end
            % in case of current density / stranded / solid conductor: overwrite
            for iEl = 1:numel(obj.CurrentPostprocessing)
                el = obj.CurrentPostprocessing(iEl).Element;
                if ismember(reg, obj.CurrentPostprocessing(iEl).Regions)
                    Jz = el.getLossDensity(t, x, y, AzPrime, reg);
                end
            end
        end

        function nuPrime = getNuPrime(obj, B)
            if obj.IsLinearMAG
                nuPrime = zeros(size(B));
            else
                B(B<=1e-8) = 1e-8; % avoid zero division
                nuPrime = (ppval(obj.HBsplineDer, B).*B - ppval(obj.HBspline, B))./(B.^2);
                % correct B values higher than measurement (assume behavior like vacuum)
                B2correct = B(B>obj.Bvalues(end));
                nuPrime(B>obj.Bvalues(end)) = (obj.Bvalues(end)-obj.Mu0*obj.Hvalues(end))/obj.Mu0./B2correct.^2; 
            end
        end

        function d2nudB2 = getd2nudB2(obj, B)
            if obj.IsLinearMAG
                d2nudB2 = zeros(size(B));
            else
                B(B<=1e-8) = 1e-8; % avoid zero division
                d2nudB2 = (ppval(obj.HBsplineDer2, B).*B.^2 - 2* ppval(obj.HBsplineDer, B).*B  + 2*ppval(obj.HBspline, B))./(B.^3);
                B2correct = B(B>obj.Bvalues(end));
                d2nudB2(B>obj.Bvalues(end)) = -2*(obj.Bvalues(end)-obj.Mu0*obj.Hvalues(end))/obj.Mu0./B2correct.^3; 
            end
        end

        function setRelativePermeability(obj, mur)
            obj.MuRel = mur;
        end

        function setConductivity(obj, sigma)
            obj.Sigma = sigma;
        end

        function fitHBspline(obj, extendBH)
            arguments
                obj 
                extendBH (1,1) logical = true
            end
            % extend B and H manually by three linear points for smooth extrapolation of vacuum permeability
            if extendBH
                deltaH = obj.Hvalues(end) - obj.Hvalues(1);
                deltaB = obj.Mu0*deltaH;
                obj.Hvalues = [obj.Hvalues, obj.Hvalues(end) + deltaH*(1:3)];
                obj.Bvalues = [obj.Bvalues, obj.Bvalues(end) + deltaB*(1:3)];
            end
            obj.HBspline = pchip(obj.Bvalues, obj.Hvalues);
            obj.HBsplineDer = fnder(obj.HBspline);
            obj.HBsplineDer2 = fnder(obj.HBsplineDer);
            obj.HBsplineInt = fnint(obj.HBspline);
        end

        function plotCharacteristics(obj)
            figure("Name", "Material Characteristics");
            % B-H curve
            subplot(2,2,1)
            obj.plotBH();
            % B-mu-curve
            subplot(2,2,2)
            obj.plotBMu();
            % B-nu-curve
            subplot(2,2,3)
            obj.plotBNu();
            % B-nuprime-curve
            subplot(2,2,4)
            obj.plotBNuPrime();
        end

        function plotBH(obj)
            b = 0:0.001:max(obj.Bvalues)*1.1;
            h = b./obj.getPermeability(b);
            % B-H curve
            plot(h, b, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.Hvalues, obj.Bvalues, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("H (A/m)")
            ylabel("B (T)")
            title("B-H-curve")
        end

        function plotBMu(obj)
            b = 0:0.001:max(obj.Bvalues)*1.1;
            mur = obj.getPermeability(b)/obj.Mu0;
            plot(b, mur, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.Bvalues, obj.Bvalues./obj.Hvalues/obj.Mu0, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Mu_r")
            title("B-mu-curve")
        end

        function plotBNu(obj)
            b = 0:0.001:max(obj.Bvalues)*1.1;
            nu = obj.getReluctivity(b);
            plot(b, nu, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.Bvalues, obj.Hvalues./obj.Bvalues, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Nu ")
            title("B-nu-curve")
        end

        function plotBNuPrime(obj)
            b = 0:0.001:max(obj.Bvalues)*1.1;
            nuP = obj.getNuPrime(b);
            plot(b, nuP, "Color", "black", LineWidth=1.5)
            grid on
            xlabel("B (T)")
            ylabel("dNu/dB ")
            title("B-nu'-curve")
        end
    end
end