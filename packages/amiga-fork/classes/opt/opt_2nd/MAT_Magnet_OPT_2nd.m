classdef MAT_Magnet_OPT_2nd < MAT_Magnet_OPT & OPT_Element_2nd
    properties (Access = public)

    end

    methods (Access = public)
        function obj = MAT_Magnet_OPT_2nd(region, regionNumbers, angle, args)
            arguments
                region 
                regionNumbers 
                angle
                args.AngleFunction = [];        % optionally if angle changes: function for the angle
                args.AngleParameter = [];       % optionally if angle changes: parameter name that is used in AngleFunction
                args.Remanence = 1;
                args.RelativePermeability = 1.05;
            end
            obj@MAT_Magnet_OPT(region, regionNumbers, angle, "AngleFunction",...
                args.AngleFunction, "AngleParameter", args.AngleParameter, "Remanence", args.Remanence, "RelativePermeability", args.RelativePermeability);
        end

        function Lxx = function_Lxx(obj, t, x, y, z)
            Lxx = zeros(obj.NumberOptimizationParameters, obj.NumberOptimizationParameters, obj.NumberOptimizationResponses);
            
            indArea = findStructName(obj.OptimizationResponses, "Area");
            if ~isempty(indArea) && ~isempty(obj.Region.dCdP)
                Lxx(:,:, indArea) = int_D_DP_DP_Omega(obj.SpaceGeo, obj.Mesh, obj.Region.dCdP(obj.SpaceGeo.gnum,:,:), obj.Region.d2CdP2(obj.SpaceGeo.gnum,:,:,:));
            end
            indVolume = findStructName(obj.OptimizationResponses, "Volume");
            if ~isempty(indVolume) && ~isempty(obj.Region.dCdP)
                Lxx(:,:, indVolume) = obj.Mesh.length*int_D_DP_DP_Omega(obj.SpaceGeo, obj.Mesh, obj.Region.dCdP(obj.SpaceGeo.gnum,:,:), obj.Region.d2CdP2(obj.SpaceGeo.gnum,:,:,:));
            end
            indCost = findStructName(obj.OptimizationResponses, "Cost");
            if ~isempty(indCost) && ~isempty(obj.Region.dCdP)
                Lxx(:,:, indCost) = obj.Mesh.length*obj.getRho()*obj.getCostperKG()*int_D_DP_DP_Omega(obj.SpaceGeo, obj.Mesh, obj.Region.dCdP(obj.SpaceGeo.gnum,:,:), obj.Region.d2CdP2(obj.SpaceGeo.gnum,:,:,:));
            end

            indTorque = findStructName(obj.OptimizationResponses, "Torque");

            if  ~isempty(obj.Region.dCdP) && ~isempty(indTorque)
                z_T = z(:, indTorque);
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);
                d2CdP2 = obj.Region.d2CdP2(obj.SpaceGeo.gnum,:,:,:);
               
                if obj.IsLinearMAG
                    nu = obj.getReluctivity();
                    dT = - obj.Mesh.length*nu*int_D_DP_DP_gradu_u_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, d2CdP2) ...
                         + obj.Mesh.length*int_D_DP_DP_gradv_Br_bot_times_z(obj.Space, obj.SpaceGeo, obj.Mesh, z_T, dCdP, d2CdP2, obj);
                else
                    error("only linear magnet implemented");
                    % dT = - obj.Mesh.length*int_D_DP_DP_gradu_u_nu_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, d2CdP2, obj, 1);
                end
                Lxx(:, :, indTorque) = dT;

                % mytimer("Br", "start")
                % change of RHS due to change of magnetization direction
                if ~isempty(obj.AngleFunction) && ~ isempty(obj.AngleParameter)
                    % Find angle index
                    indAngle = findStructName(obj.OptimizationParameters, obj.AngleParameter);
                    angle = obj.AngleFunction(x(indAngle));
                    % Numerical evaluation how parameter changes angle
                    Step = (obj.OptimizationParameters(indAngle).Max - obj.OptimizationParameters(indAngle).Min)*1e-4;
                    dAngledP = (obj.AngleFunction(x(indAngle)+Step) - angle)/Step;
                    d2AngledP2 = (obj.AngleFunction(x(indAngle)+Step) - 2*angle + obj.AngleFunction(x(indAngle)-Step))/(Step^2);
                    
                    % d_alpha_d_alpha
                    cont_1 = - obj.Mesh.length*dAngledP^2*z_T'*int_gradv_Br_bot(obj.Space, obj.Mesh, obj);
                    cont_2 =   obj.Mesh.length*d2AngledP2.*z_T'*int_d_dalpha_gradv_Br_bot(obj.Space, obj.Mesh, obj);
                    % d_alpha_dP
                    cont_3 = obj.Mesh.length*dAngledP*z_T'*int_D_DP_d_alpha_gradv_Br_bot(obj.Space, obj.SpaceGeo, obj.Mesh, dCdP, obj);

                    Lxx(indAngle, indAngle, indTorque) = Lxx(indAngle, indAngle, indTorque) + cont_1 + cont_2 + cont_3(indAngle);
                end
                % mytimer("Br", "stop")
                
            end
        end

        function Lxy = function_Lxy(obj, t, x, y, z)
            % Lxy = zeros(obj.NumberOptimizationParameters, obj.Space.ndof, obj.NumberOptimizationResponses);
            Lxy = sptensor([], [], [obj.NumberOptimizationParameters, obj.Space.ndof, obj.NumberOptimizationResponses]);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(obj.Region.dCdP) && ~isempty(indTorque)
                z_T = z(:, indTorque);
                dCdP = obj.Region.dCdP(obj.SpaceGeo.gnum,:,:);

                if obj.IsLinearMAG
                    nu = obj.getReluctivity();
                    dT = - obj.Mesh.length*nu*int_D_DP_gradu_u_gradv(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, z_T, dCdP);
                    % dT = dT - nu/obj.Length*int_D_DP_gradu_u_gradv_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                    %         obj.SpacesGeo, obj.SpacesGeoEval, obj.Meshes, obj.MeshesEval, z_T, obj.dCdP, patches)';
                else
                    % first contribution
                    dT = - obj.Mesh.length*int_D_DY_DP_gradu_u_nu_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, obj, 1) ...;
                         - obj.Mesh.length*int_D_DP_gradu_nu_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, obj, 1);
                    % dT = dT - int_D_DY_DP_gradu_u_nu_gradv_v_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                    %     obj.SpacesGeo, obj.SpacesGeoEval, obj.Meshes, obj.MeshesEval, y, z_T, obj.dCdP, mat, obj.Length, patches)';
                    % second contribution
                    % dT = dT - int_D_DP_gradu_nu_gradv_v_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                    %     obj.SpacesGeo, obj.SpacesGeoEval, obj.Meshes, obj.MeshesEval, y, z_T, obj.dCdP, mat, obj.Length, patches)';
                end
                Lxy(:, :, indTorque) = sptensor(dT');
                % Lxy(:, :, indTorque) = dT';
            end
        end

        function Lyy = function_Lyy(obj, t, x, y, z)
            if obj.NumberOptimizationResponses == 0
                Lyy = zeros(obj.Space.ndof, obj.Space.ndof, obj.NumberOptimizationResponses);
                return;
            end
            
            Lyy = sptensor([],[], [obj.Space.ndof, obj.Space.ndof, obj.NumberOptimizationResponses]);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(indTorque)
                z_T = z(:, indTorque);
                if ~obj.IsLinearMAG
                    Lyy(:, :, indTorque) = sptensor(- obj.Mesh.length*int_D_DY_DY_gradu_u_nu_gradv_v(obj.Space, obj.Space, obj.Mesh, y, z_T, obj, 1));
                end
            end
        end
        
    end
end