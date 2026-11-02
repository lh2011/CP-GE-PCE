classdef MAT_Thermal < MAT & EQS_Element
    properties
        IsLinearTH = true;
        HeatCapacity % Volumetric [W/(kg*K)]
        Conductivity % [J/(m^3*K)]
        Tvalues, Kvalues;
        TKspline, TKsplineDer, TKsplineInt;
    end

    methods
        function obj = MAT_Thermal(region, regionNumbers, args)
            arguments
                region;
                regionNumbers;
                args.FEMMname = [];
                args.Conductivity = [];
                args.HeatCapacity = [];
                args.PlotColor = [0.9,0.9,0.9];
            end
            obj.Region = region;
            obj.RegionNumbers = regionNumbers;
            
            obj.PlotColor = args.PlotColor;
            % condensing MESH and SPACE coming from the geometry class
            elementNumbers = region.getPatchNumbers(regionNumbers);
            obj.Mesh = MESH_2D(region.MeshesEval(elementNumbers), "length", region.Length);
            obj.Space = SPACE_2D(region.SpacesEval(elementNumbers), region.Spaces.gnum(elementNumbers), region.Spaces.ncomp);
            obj.NumberDOF = 0;                  % material has no intrinsic dofs but writes to the external space object
            region.setMaterial(obj, regionNumbers);    % assign material to region

            obj.setConductivity(args.Conductivity);
            obj.setHeatCapacity(args.HeatCapacity);
            if ~isempty(args.FEMMname)
                obj.loadFEMMmaterial(args.FEMMname);
            end
        end

        function cp = getHeatCapacity(obj)
            cp = obj.HeatCapacity;
        end

        function cp = setHeatCapacity(obj, cp)
            obj.HeatCapacity = cp;
        end
        
        function kappa = getConductivity(obj)
            kappa = obj.Conductivity;
        end

        function setLinear(obj, trueFalse)
            obj.IsLinearTH = trueFalse;
        end

        function kappa = setConductivity(obj, kappa)
            obj.Conductivity = kappa;
        end

        function postprocess(obj, t, y, yprime)
            obj.Solution.Time = t;
            obj.Solution.Temperature = y;
            obj.Solution.TemperatureDer = yprime;
        end

        function loadFEMMmaterial(obj, femmName)
            obj.Name = femmName;
            materialStruct = parseFEMMthermal(femmName);
            if isempty(materialStruct)
                error('MAT_Thermal: FEMM Material not found');
            end
            obj.setConductivity(materialStruct.Kx);
            obj.setHeatCapacity(materialStruct.Kt*1e6);
            if isempty(materialStruct.TK)
                obj.IsLinearTH = true;
            else
                obj.IsLinearTH = false;
                obj.Tvalues = materialStruct.TK(:,1)';
                obj.Kvalues = materialStruct.TK(:,2)';
                obj.fitTKspline();
            end
        end

        function fitTKspline(obj)
            obj.TKspline = pchip(obj.Tvalues, obj.Kvalues);
            obj.TKsplineDer = fnder(obj.TKspline);
            obj.TKsplineInt = fnint(obj.TKspline);
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            if nargout == 2 
                Mlin = obj.Mesh.length*int_u_v(obj.Space, obj.Space, obj.Mesh, obj.getHeatCapacity());
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.Space.ndof, obj.Space.ndof);
            if obj.IsLinearTH
                if nargout == 2
                    Klin = - obj.Mesh.length*int_gradu_gradv(obj.Space, obj.Space, obj.Mesh, obj.getConductivity());
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
            if obj.IsLinearTH
                J = sparse(obj.Space.ndof, obj.Space.ndof);
            else
                error("Nonlinear material not yet implemented")
            end
        end
    end
end