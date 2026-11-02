% class holding the meshes and spaces of a boundary of one region
classdef IGA_CouplingInterface < handle
    properties (Access = public)
        Region;                         % Region reference
        PatchBoundaries;                % Boundary numbers of region
        HarmonicsSin, HarmonicsCos;     % Used sin and cos values for coupling 
        HarmonicsAll;                   % Union of sin/cos for efficient rotation evaluation
        CouplingIndices;                % Indices for used sin/cos harmonics in all harmonics
        Radius;                         % Radius of Coupling Interface
        ArcLength;                      % Interface ArcLength
        Spaces, SpacesEval;             % 1D boundary spaces
        MeshesEval;                     % 1D boundary meshes
        NumberBoundaryDOFs;             % Number of FEM dofs that are coupled

        CouplingMatrixFull;             % Unrotated Matrix of all sin/cos/sin/cos... harmonics
        CouplingMatrixInit;             % Unrotated Matrix of used sin/cos harmonics
        CouplingMatrix;                 % Rotated   Matrix of used sin/cos harmonics
        CouplingMatrixDer;              % Derivative of CouplingMatrix w.r.t rotation angle
        CouplingMatrixDer2;             % Second derivative of CouplingMatrix w.r.t rotation angle

        RotationAngle = 0;              % Angle shift for evaluation
        RotationMatrix;                 % Matrix to gain CouplingMatrix from CouplingMatrixInit
        RotationMatrixDer;              % Derivative of RotationMatrix w.r.t rotation angle
        RotationMatrixDer2;             % Second derivative of RotationMatrix w.r.t rotation angle
    end

    methods (Access = public)
        function obj = IGA_CouplingInterface(region, patchBoundaries, nGaussPoints)
            obj.Region = region;
            obj.PatchBoundaries = patchBoundaries;
            obj.NumberBoundaryDOFs = numel(region.getBoundaryDOFs(patchBoundaries));
            if isempty(nGaussPoints)
                % By default, take much more quadrature points for better sin/cos integral evaluation
                nGaussPoints = region.DegreeQuadrature*10;
            end
            obj.generateSpacesAndMeshes(nGaussPoints);
        end

        function generateSpacesAndMeshes(obj, nGaussNodes)
            % region1
            for iBnd = 1:numel(obj.PatchBoundaries)
                patch = obj.Region.PatchBoundaries(obj.PatchBoundaries(iBnd)).patches;
                side = obj.Region.PatchBoundaries(obj.PatchBoundaries(iBnd)).faces;
                if ~exist("nGaussNodes", "var")
                    % use quadrature rule of region:
                    obj.MeshesEval{iBnd} = msh_eval_boundary_side(obj.Region.Meshes.msh_patch{patch}, side);
                    msh_side_int = msh_boundary_side_from_interior(obj.Region.Meshes.msh_patch{patch}, side);
                else
                    % Recalculate if different quadrature is wanted on interface
                    [~, zeta] = kntrefine (obj.Region.Geometry(patch).nurbs.knots, obj.Region.SubdivisionsPatches-1, obj.Region.FormFunctionDegree, obj.Region.FormFunctionDegree - 1);
                    rule      = msh_gauss_nodes (nGaussNodes);
                    [qn, qw]  = msh_set_quad_nodes (zeta, rule);
                    mesh = msh_cartesian (zeta, qn, qw, obj.Region.Geometry(patch));
                    obj.MeshesEval{iBnd} = msh_eval_boundary_side(mesh, side);
                    msh_side_int = msh_boundary_side_from_interior(mesh, side);
                end
                obj.Spaces{iBnd} = obj.Region.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval{iBnd} = sp_precompute(obj.Spaces{iBnd}, msh_side_int, 'value', true, 'gradient', true);
            end

            obj.ArcLength = 0;
            for iBnd = 1:numel(obj.PatchBoundaries)
                obj.ArcLength = obj.ArcLength + int_Omega(obj.MeshesEval{iBnd});
            end

            obj.Radius = squeeze((obj.MeshesEval{1}.geo_map(1,1,1).^2 + obj.MeshesEval{1}.geo_map(2,1,1).^2).^0.5);

            % check that coupling boundary is a circle
            for iBnd = 1:numel(obj.PatchBoundaries)
                meshRadius = squeeze((obj.MeshesEval{iBnd}.geo_map(1,:,:).^2 + obj.MeshesEval{iBnd}.geo_map(2,:,:).^2).^0.5);
                assert(all(abs(meshRadius-obj.Radius) < 1e-10, "all"), "Coupling Boundary must be a circle!");
            end
        end

        function setHarmonics(obj, sinValues, cosValues)
            obj.HarmonicsSin = reshape(unique(sinValues), 1, []);
            obj.HarmonicsCos = reshape(unique(cosValues), 1, []);
            obj.HarmonicsAll = union(obj.HarmonicsSin, obj.HarmonicsCos);

            [~, indexSin, ~] = intersect(obj.HarmonicsAll, obj.HarmonicsSin);
            [~, indexCos, ~] = intersect(obj.HarmonicsAll, obj.HarmonicsCos);
            obj.CouplingIndices = union(2*indexSin-1, 2*indexCos);
            obj.generateCouplingMatrices();
        end
        
        function generateCouplingMatrices(obj)
            % define array functions for all harmonics
            SinPhi =  @(x, y) sin(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));
            CosPhi =  @(x, y) cos(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));

            obj.CouplingMatrixFull = spalloc (obj.Region.NumberDOF, 2*numel(obj.HarmonicsAll), 2*numel(obj.HarmonicsAll)*numel(obj.Region.getBoundaryDOFs(obj.PatchBoundaries)));
            for iBnd = 1:numel(obj.PatchBoundaries)
                patchDOFs = obj.Region.Spaces.gnum{obj.Region.PatchBoundaries(obj.PatchBoundaries(iBnd)).patches};
                obj.CouplingMatrixFull(patchDOFs, 1:2:end) = obj.CouplingMatrixFull (patchDOFs, 1:2:end) + int_fs_v (obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, SinPhi);
                obj.CouplingMatrixFull(patchDOFs, 2:2:end) = obj.CouplingMatrixFull (patchDOFs, 2:2:end) + int_fs_v (obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, CosPhi);
            end
            obj.CouplingMatrixInit = obj.CouplingMatrixFull(:, obj.CouplingIndices);
            obj.setRotationAngle(obj.RotationAngle);
        end

        function setRotationAngle(obj, angleRad)
            obj.RotationAngle = angleRad;

            nHarm = numel(obj.HarmonicsAll);
            rows = [1:2:nHarm*2, 1:2:nHarm*2, 2:2:nHarm*2, 2:2:nHarm*2];
            cols = [1:2:nHarm*2, 2:2:nHarm*2, 1:2:nHarm*2, 2:2:nHarm*2];
            vals = [cos(obj.HarmonicsAll.*obj.RotationAngle), sin(obj.HarmonicsAll.*obj.RotationAngle), ...
                   -sin(obj.HarmonicsAll.*obj.RotationAngle), cos(obj.HarmonicsAll.*obj.RotationAngle)];
            R = sparse(rows, cols, vals, 2*nHarm, 2*nHarm);

            valsDer = [-obj.HarmonicsAll.*sin(obj.HarmonicsAll.*obj.RotationAngle), obj.HarmonicsAll.*cos(obj.HarmonicsAll.*obj.RotationAngle), ...
                       -obj.HarmonicsAll.*cos(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.*sin(obj.HarmonicsAll.*obj.RotationAngle)];
            Rder = sparse(rows, cols, valsDer, 2*nHarm, 2*nHarm);

            valsDer2 = [-obj.HarmonicsAll.^2.*cos(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.^2.*sin(obj.HarmonicsAll.*obj.RotationAngle), ...
                         obj.HarmonicsAll.^2.*sin(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.^2.*cos(obj.HarmonicsAll.*obj.RotationAngle)];
            Rder2 = sparse(rows, cols, valsDer2, 2*nHarm, 2*nHarm);

            obj.RotationMatrix = R(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer = Rder(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer2 = Rder2(obj.CouplingIndices, obj.CouplingIndices);

            obj.CouplingMatrix = obj.CouplingMatrixFull*R(:, obj.CouplingIndices);
            obj.CouplingMatrixDer = obj.CouplingMatrixFull*Rder(:, obj.CouplingIndices);
            obj.CouplingMatrixDer2 = obj.CouplingMatrixFull*Rder2(:, obj.CouplingIndices);
        end

        function plotCouplingInterface(obj, col)
            if ~exist("color", "var")
                col = TUDa_getColor("9b");
            end
            lw = 2;
            for iBnd = obj.PatchBoundaries
                plot(obj.Region.PatchBoundaries(iBnd).PlotPoints{1}, obj.Region.PatchBoundaries(iBnd).PlotPoints{2}, "Color", col, "LineWidth", lw);
            end
        end
    end
end