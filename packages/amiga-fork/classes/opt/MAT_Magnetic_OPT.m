classdef MAT_Magnetic_OPT < MAT_Magnetic & OPT_Element
    properties
        SpaceGeo
    end

    methods (Access = public)
        function obj = MAT_Magnetic_OPT(region, regionNumbers, args)
            arguments
                region;
                regionNumbers;
                args.FEMMname = [];
                args.LaminationThickness = 0;
                args.FillFactor = 1;
                args.PlotColor = [0.9,0.3,0.4];
            end
            args = namedargs2cell(args);
            obj@MAT_Magnetic(region, regionNumbers, args{:});
        end

        function optimizationUpdate(obj, x)
            % update spaces and meshes for new geometry
            elementNumbers = obj.Region.getPatchNumbers(obj.RegionNumbers);
            obj.Mesh = MESH_2D(obj.Region.MeshesEval(elementNumbers), "length", obj.Region.Length);
            obj.Space = SPACE_2D(obj.Region.SpacesEval(elementNumbers), obj.Region.Spaces.gnum(elementNumbers), obj.Region.Spaces.ncomp);
            obj.SpaceGeo = SPACE_2D(obj.Region.SpacesGeoEval(elementNumbers), obj.Region.SpacesGeo.gnum(elementNumbers), obj.Region.Spaces.ncomp);
            
            for iExc = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iExc).Excitation, "OPT_Element")
                    obj.Excitations(iExc).Excitation.optimizationUpdate(x);
                end
            end
        end

        function r = responses(obj, t, x, y)
            r = zeros(obj.NumberOptimizationResponses, 1);

            indCost = findStructName(obj.OptimizationResponses, "Cost");
            if ~isempty(indCost)
                r(indCost) = obj.getCost();
            end
            indVolume = findStructName(obj.OptimizationResponses, "Volume");
            if ~isempty(indVolume)
                r(indVolume) = obj.getVolume();
            end
            indArea = findStructName(obj.OptimizationResponses, "Area");
            if ~isempty(indArea)
                r(indArea) = obj.getArea();
            end
        end

        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = sparse(obj.NumberOptimizationResponses, obj.Space.ndof);
        end

        function drdx = partial_r_partial_x(obj, t, x, y)
            drdx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            
            indCost = findStructName(obj.OptimizationResponses, "Cost");
            if ~isempty(indCost) && ~isempty(obj.Region.dCdP)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                drdx(indCost, :) = obj.getCostperKG()*obj.getRho()*obj.Mesh.length*int_D_DP_Omega(obj.SpaceGeo, obj.Mesh, dCdP)';
            end
            indVolume = findStructName(obj.OptimizationResponses, "Volume");
            if ~isempty(indVolume) && ~isempty(obj.Region.dCdP)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                drdx(indVolume, :) = obj.Mesh.length*int_D_DP_Omega(obj.SpaceGeo, obj.Mesh, dCdP)';
            end
            indArea = findStructName(obj.OptimizationResponses, "Area");
            if ~isempty(indArea) && ~isempty(obj.Region.dCdP)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                drdx(indArea, :) = int_D_DP_Omega(obj.SpaceGeo, obj.Mesh, dCdP)';
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            if ~isempty(obj.Region.dCdP)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                if obj.IsLinearMAG == true
                    nu = obj.getReluctivity();
                    dedx = - obj.Mesh.length*nu*int_D_DP_gradu_u_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, dCdP);
                else
                    dedx = - obj.Mesh.length*int_D_DP_gradu_u_nu_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, dCdP, obj);
                end
            else
                dedx = zeros(obj.Space.ndof, obj.NumberOptimizationParameters);
            end
            
            for iEXC = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iEXC).Excitation, "OPT_Element")
                    dedx = dedx + obj.Excitations(iEXC).Excitation.partial_e_partial_x(t, x, y);
                end
            end
        end

        function Lx = function_Lx(obj, t, x, y, z)
            Lx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(obj.Region.dCdP) && ~isempty(indTorque)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                z_Torque = z(:,indTorque);
                if obj.IsLinearMAG == true
                    nu = obj.getReluctivity();
                    Lx(indTorque, :) = - obj.Mesh.length*nu*int_D_DP_z_grady_y_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_Torque, dCdP);
                else
                    Lx(indTorque, :) = - obj.Mesh.length*int_D_DP_z_grady_y_nu_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_Torque, dCdP, obj);
                end
            end

            for iExc = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iExc).Excitation, "OPT_Element")
                    Lx = Lx + obj.Excitations(iExc).Excitation.function_Lx(t, x, y, z);
                end
            end
        end
    end
end