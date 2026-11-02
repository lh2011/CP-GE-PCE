classdef MAT_Mechanical < MAT & EQS_Element
    properties
        IsLinearMECH = true;
        % StressValues, StrainValues;
        % SEspline, SEsplineDer, SEsplineInt;
        PoissonRatio, YoungsModulus;
    end

    methods
        function obj = MAT_Mechanical(region, regionNumbers, args)
            arguments
                region;
                regionNumbers;
                args.PoissonRatio = [];
                args.ElasticityModulus = [];
                args.PlotColor = [0.9,0.9,0.9];
            end
            obj.Region = region;
            obj.RegionNumbers = regionNumbers;
            
            obj.PlotColor = args.PlotColor;
            % condensing MESH and SPACE coming from the geometry class
            elementNumbers = region.getPatchNumbers(regionNumbers);
            obj.Mesh = MESH_2D(region.MeshesEval(elementNumbers), "length", region.Length);
            obj.Space = SPACE_2D(region.SpacesEval(elementNumbers), region.Spaces.gnum(elementNumbers), region.Spaces.ncomp);
            obj.NumberDOF = 0;                         % material has no intrinsic dofs but writes to the external space object
            region.setMaterial(obj, regionNumbers);    % assign material to region

            obj.setPoissonRatio(args.PoissonRatio);
            obj.setYoungsModulus(args.ElasticityModulus);
        end

        function E = getYoungsModulus(obj)
            E = obj.YoungsModulus;
        end

        function E = setYoungsModulus(obj, E)
            obj.YoungsModulus = E;
        end
        
        function nu = getPoissonRatio(obj)
            nu = obj.PoissonRatio;
        end

        function nu = setPoissonRatio(obj, nu)
            obj.PoissonRatio = nu;
        end

        function setLinear(obj, trueFalse)
            obj.IsLinearMECH = trueFalse;
        end     

        function stress = getStress(obj, strain)
            arguments
                obj 
                strain (3,:) double;
            end
            % lets just assume plain stress for now
            nu = obj.getPoissonRatio();
            E = obj.getYoungsModulus();
            C = E/(1-nu^2)* [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
            stress = tensorprod(C, strain, 2, 1);
        end

        function postprocess(obj, t, y, yprime)
            obj.Solution.Time = t;
        end

        function fitSEspline(obj)
            obj.SEspline = pchip(obj.StressValues, obj.StrainValues);
            obj.SEsplineDer = fnder(obj.SEspline);
            obj.SEsplineInt = fnint(obj.SEspline);
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            Mlin = sparse(obj.Space.ndof, obj.Space.ndof);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            if obj.IsLinearMECH
                if nargout == 2
                    nu = obj.getPoissonRatio();
                    E = obj.getYoungsModulus();
                    % plain stress
                    C = E/(1-nu^2)* [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
                    Klin = - obj.Mesh.length*int_eu_C_ev(obj.Space, obj.Space, obj.Mesh, C);
                end
            else
                error("nonlinear material not yet implemented")
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = zeros(obj.Space.ndof, 1);
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
            if obj.IsLinearMECH
                J = sparse(obj.Space.ndof, obj.Space.ndof);
            else
                error("Nonlinear material not yet implemented")
            end
        end
    end
end