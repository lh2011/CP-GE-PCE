classdef MAT_Magnet_OPT < MAT_Magnet & OPT_Element
    properties (Access = public)
        SpaceGeo
        AngleFunction;
        AngleParameter;
    end

    methods (Access = public)
        function obj = MAT_Magnet_OPT(region, regionNumbers, angle, args)
            arguments
                region 
                regionNumbers 
                angle
                args.AngleFunction = [];        % optionally if angle changes: function for the angle
                args.AngleParameter = [];       % optionally if angle changes: parameter name that is used in AngleFunction
                args.Remanence = 1;
                args.RelativePermeability = 1.05;
            end
            obj@MAT_Magnet(region, regionNumbers, angle, "Remanence", args.Remanence, "RelativePermeability", args.RelativePermeability);
            obj.AngleFunction = args.AngleFunction;
            obj.AngleParameter = args.AngleParameter;
        end

        function optimizationUpdate(obj, x)
            % update spaces and meshes for new geometry
            elementNumbers = obj.Region.getPatchNumbers(obj.RegionNumbers);
            obj.Mesh = MESH_2D(obj.Region.MeshesEval(elementNumbers), "length", obj.Region.Length);
            obj.Space = SPACE_2D(obj.Region.SpacesEval(elementNumbers), obj.Region.Spaces.gnum(elementNumbers), obj.Region.Spaces.ncomp);
            obj.SpaceGeo = SPACE_2D(obj.Region.SpacesGeoEval(elementNumbers), obj.Region.SpacesGeo.gnum(elementNumbers), obj.Region.Spaces.ncomp);
            % update remanence direction
            if ~isempty(obj.AngleFunction) && ~isempty(obj.AngleParameter)
                indAngle = findStructName(obj.OptimizationParameters, obj.AngleParameter);
                angle = obj.AngleFunction(x(indAngle));
                obj.Angle = @(x,y) ones(size(x))* angle;
            end
        end

        function r = responses(obj, t, x, y)
            r = zeros(obj.NumberOptimizationResponses, 1);
            if obj.NumberOptimizationResponses == 0 
                return 
            end
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

        function Lx = function_Lx(obj, t, x, y, z)
            Lx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(obj.Region.dCdP) && ~isempty(indTorque)
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum, :, :);
                z_Torque = z(:, indTorque);
                nu = obj.getReluctivity();
                Lx(indTorque, :) = - obj.Mesh.length*nu*int_D_DP_z_grady_y_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_Torque, dCdP)...
                                   + obj.Mesh.length*int_D_DP_z_gradv_Br_bot(obj.Space, obj.SpaceGeo, obj.Mesh, z_Torque, dCdP, obj);
                
                % change of RHS due to change of magnetization direction
                if ~isempty(obj.AngleFunction) && ~isempty(obj.AngleParameter)
                    indAngle = findStructName(obj.OptimizationParameters, obj.AngleParameter);
                    angle = obj.AngleFunction(x(indAngle));

                    Step = (obj.OptimizationParameters(indAngle).Max - obj.OptimizationParameters(indAngle).Min)*1e-4;
                    dAngledP = (obj.AngleFunction(x(indAngle)+Step) - angle)/Step;

                    Lx(indTorque, indAngle) = Lx(indTorque, indAngle) + obj.Mesh.length* dAngledP.* z_Torque'* int_d_dalpha_gradv_Br_bot(obj.Space, obj.Mesh, obj);
                end
            end

            for iExc = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iExc).Excitation, "OPT_Element")
                    Lx = Lx + obj.Excitations(iExc).Excitation.function_Lx(t, x, y, z);
                end
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            if isempty(obj.Region.dCdP)
                dedx = zeros(obj.Space.ndof, obj.NumberOptimizationParameters);
                return;
            end
            dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
            if obj.IsLinearMAG == true
                nu = obj.getReluctivity();
                dedx = - obj.Mesh.length*nu*int_D_DP_gradu_u_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, dCdP);
            else
                dedx = - obj.Mesh.length*int_D_DP_gradu_u_nu_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, dCdP, obj);
            end

            dedx = dedx + obj.Mesh.length*int_D_DP_gradv_Br_bot(obj.Space, obj.SpaceGeo, obj.Mesh, dCdP, obj);
            
            % change of RHS due to change of magnetization direction
            if ~isempty(obj.AngleFunction) && ~ isempty(obj.AngleParameter)
                % Find angle index
                indAngle = findStructName(obj.OptimizationParameters, obj.AngleParameter);
                angle = obj.AngleFunction(x(indAngle));
                % Numerical evaluation how parameter changes angle
                Step = (obj.OptimizationParameters(indAngle).Max - obj.OptimizationParameters(indAngle).Min)*1e-4;
                dAngledP = (obj.AngleFunction(x(indAngle)+Step) - angle)/Step;
                dedx(:, indAngle) = dedx(:, indAngle) + obj.Mesh.length*dAngledP*int_d_dalpha_gradv_Br_bot(obj.Space, obj.Mesh, obj);
            end
        end
    end
end