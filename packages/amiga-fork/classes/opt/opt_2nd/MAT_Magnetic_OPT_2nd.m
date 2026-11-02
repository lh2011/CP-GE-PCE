classdef MAT_Magnetic_OPT_2nd < MAT_Magnetic_OPT & OPT_Element_2nd
    
    methods (Access = public)
        function obj = MAT_Magnetic_OPT_2nd(region, regionNumbers, args)
            arguments
                region;
                regionNumbers;
                args.FEMMname = [];
                args.LaminationThickness = 0;
                args.FillFactor = 1;
                args.PlotColor = [0.9,0.3,0.4];
            end
            args = namedargs2cell(args);
            obj@MAT_Magnetic_OPT(region, regionNumbers, args{:});
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
                    dT = - obj.Mesh.length*nu*int_D_DP_DP_gradu_u_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, d2CdP2);
                else
                    dT = - obj.Mesh.length*int_D_DP_DP_gradu_u_nu_gradv_v(obj.Space, obj.Space, obj.SpaceGeo, obj.Mesh, y, z_T, dCdP, d2CdP2, obj);
                end
                
                Lxx(:, :, indTorque) = dT;
            end

            for iExc = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iExc).Excitation, "OPT_Element_2nd")
                    Lxx = Lxx + obj.Excitations(iExc).Excitation.function_Lxx(t, x, y, z);
                end
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