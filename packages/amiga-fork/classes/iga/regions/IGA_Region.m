classdef IGA_Region < EQS_System

    properties
        % quality of computation
        FormFunctionDegree ;                % degree of form functions 1-linear, 2-quadratic, 3-cubic, ...
        DegreeQuadrature;                   % degree of quadrature, usually use FormFunctionDegree
        SubdivisionsPatches;                % subdivision for each patch

        % domain definitions
        Length;                             % axial length
        Surface;                            % NURBS representation of single patches
        Geometry;                           % parametrization of the mapping to the physical domain and derivatives
        PatchBoundaries;                         % information about the patches with an outer boundary
        PatchInterfaces;                         % information about the connections between patches
        BoundaryInterfaces;                 % boundary information
        ControlPoints;                      % control points and weights
        NumberPatches;                      % patch counter

        % computation of spaces and meshes
        Spaces, SpacesEval;                 % space stores the shape functions of discrete spaces
        SpacesGeo, SpacesGeoEval;           % space of mapping functions for the geometry
        Meshes, MeshesEval;                 % mesh stores partition of domain and quadrature rule

        NumberRegions, NumberBoundaries;
        Regions = struct('Name', {}, 'Patches', {}, 'Material', {});
        Boundaries = struct('Name', {}, 'Patches', {}, 'Sides', {});
        
        EvaluationPointsParam, EvaluationPatches;       % parametric points and values for postprocessing evaluation
        EvaluationShapeValues;        % shape values
        PlotValuesNeedUpdate = true;

        % plotting
        PlotResolution = [5, 5];
        RotationAngle = 0;      % current rotation angle of geometry
        PlotAngle = 0;          % angle when geometry was last plotted, if ~=RotationAngle, update plot points
        PlotTool;
    end

    methods (Abstract)
        generateSpaces(obj)
    end

    methods
        function obj = IGA_Region(subdivisions, degrees, length, args)
            arguments
                subdivisions (1, 2) double  = [5, 5];
                degrees (1, 2) double = [2, 2];
                length (1, 1) double = 1;
                args.DegreeQuadrature = degrees + [1, 1];
            end
            obj.SubdivisionsPatches = subdivisions;
            obj.FormFunctionDegree = degrees;
            obj.Length = length;
            obj.DegreeQuadrature = args.DegreeQuadrature;

            obj.resetMaterials();
            obj.resetBoundaryConditions();

            obj.PlotTool = PLT_RectangularPattern();
        end

        function importSurface(obj, surface, regions, boundaries, args)
            arguments
                obj 
                surface 
                regions = [];
                boundaries = [];
                args.Tolerance = 1e-8;
            end
            if isempty(obj.Name)
                obj.Name = inputname(1);
            end
            obj.Surface = surface;
            warning('off', 'nrbderiv:SecondDerivative')
            warning('off', 'geopdes:nrbmultipatch')
            [obj.Geometry, obj.PatchBoundaries, obj.PatchInterfaces, ~, obj.BoundaryInterfaces] = mp_geo_load(obj.Surface, args.Tolerance);
            obj.NumberPatches = numel(obj.Geometry);
            % if no regions are specified: each patch is own region
            if isempty(regions)
                for iPatch = 1:obj.NumberPatches
                    regions.("Patch"+iPatch) = iPatch;
                end
            end

            regionNames = fieldnames(regions);
            obj.NumberRegions = numel(regionNames);
            % add region information
            for iRegion = 1:numel(regionNames)
                obj.Regions(iRegion).Name = regionNames{iRegion};
                obj.Regions(iRegion).Patches = regions.(regionNames{iRegion});
                for patch = obj.Regions(iRegion).Patches 
                    obj.Geometry(patch).Region = iRegion;
                end
            end
            % check material interfaces
            for iInt = 1:numel(obj.PatchInterfaces)
                if obj.Geometry(obj.PatchInterfaces(iInt).patch1).Region ~= obj.Geometry(obj.PatchInterfaces(iInt).patch2).Region
                    obj.PatchInterfaces(iInt).IsMaterialInterface = true;
                else
                    obj.PatchInterfaces(iInt).IsMaterialInterface = false;
                end
            end

            % if no boundaries are specified: each patch boundary is own boundary
            if isempty(boundaries)
                for iBnd = 1:numel(obj.PatchBoundaries)
                    boundaries.("Boundary"+iBnd).patches = obj.PatchBoundaries(iBnd).patches;
                    boundaries.("Boundary"+iBnd).sides = obj.PatchBoundaries(iBnd).faces;
                end
            end
            boundaryNames = fieldnames(boundaries);
            obj.NumberBoundaries = numel(boundaryNames);
            % add region information
            for iBnd = 1:numel(boundaryNames)
                obj.Boundaries(iBnd).Name = boundaryNames{iBnd};
                obj.Boundaries(iBnd).Patches = boundaries.(boundaryNames{iBnd}).patches;
                obj.Boundaries(iBnd).Sides = boundaries.(boundaryNames{iBnd}).sides;
                obj.Boundaries(iBnd).PatchBoundary = obj.getBoundaryNumbers(obj.Boundaries(iBnd).Patches, obj.Boundaries(iBnd).Sides);
            end
            mytimer("generateMeshes", "start");
            obj.generateMeshes();
            mytimer("generateMeshes", "stop");
            mytimer("generateSpaces", "start");
            obj.generateSpaces();
            mytimer("generateSpaces", "stop");

            obj.NumberDOF = obj.Spaces.ndof;
        end

        function generateMeshes(obj)
            progressbar([obj.Name ': Generating meshes: '])
            % pre-evaluated spaces and meshes
            % obj.MeshesEval = struct();%[];%cell (1, obj.NumberPatches);
            obj.MeshesEval = cell (1, obj.NumberPatches);
            obj.SpacesGeoEval = cell (1, obj.NumberPatches);
            % temporary cell arrays
            meshes = cell (1, obj.NumberPatches);
            spacesGeo = cell (1, obj.NumberPatches);

            % create geometry information
            for iPatch = 1:obj.NumberPatches
                % physical mapping and quadrature
                [~, zeta] = kntrefine (obj.Geometry(iPatch).nurbs.knots, obj.SubdivisionsPatches-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
                rule      = msh_gauss_nodes (obj.DegreeQuadrature);
                [qn, qw]  = msh_set_quad_nodes (zeta, rule);
                meshes{iPatch} = msh_cartesian (zeta, qn, qw, obj.Geometry(iPatch));
                spacesGeo{iPatch} = sp_nurbs(obj.Geometry(iPatch).nurbs, meshes{iPatch});

                % mesh and mapping evaluation
                obj.MeshesEval{iPatch} = msh_precompute(meshes{iPatch});
                obj.MeshesEval{iPatch}.length = obj.Length;
                
                obj.SpacesGeoEval{iPatch}  = sp_precompute(spacesGeo{iPatch}, meshes{iPatch}, 'value', true, 'gradient', true);

                progressbar(iPatch/obj.NumberPatches);
            end

            obj.Meshes = msh_multipatch(meshes, obj.PatchBoundaries);
            obj.SpacesGeo = sp_multipatch(spacesGeo, obj.Meshes, obj.PatchInterfaces, obj.BoundaryInterfaces);

            %  update control points
            obj.ControlPoints = zeros(obj.SpacesGeo.ndof, 4);
            globalNumbering = obj.SpacesGeo.gnum;   % fetch once to save time. Matlab does some strange copying otherwise
            for iPatch = 1:obj.NumberPatches
                localIndex = globalNumbering{iPatch};
                % localIndex = obj.SpacesGeo.gnum{iPatch};  % instead of this, it takes long time for some reason
                obj.ControlPoints(localIndex, 1) = reshape(obj.Geometry(iPatch).nurbs.coefs(1, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 2) = reshape(obj.Geometry(iPatch).nurbs.coefs(2, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 3) = reshape(obj.Geometry(iPatch).nurbs.coefs(3, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 4) = reshape(obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
            end
            obj.calculateArea();
            obj.calculatePlotPoints();
            obj.calculatePlotPointsVector();
            progressbar('finished');
        end

        function resetMaterials(obj)
            obj.Klin = []; obj.Mlin = []; obj.Flin = [];
            for iRegion = 1:obj.NumberRegions
                obj.Regions(iRegion).Material = [];
            end
            obj.resetElements();
        end

        function setMaterial(obj, material, regions)
            arguments
                obj 
                material 
                regions (1,:) double;
            end
            if ~ismember(regions, 1:obj.NumberRegions)
                error("IGA_Region:NotValidRegion: Region number %i is not with %i number of regions", regions, obj.NumberRegions);
            end
            % if ~isempty(obj.Regions(regions).Material)
            %     warning("IGA_Region:DuplicateMaterial: Material for region %i was overwritten", regions);
            %     "TBD really overwrite material, remove element or overwrite element"
            %     % resetMaterials(region)
            % end
            for reg = regions
                obj.Regions(reg).Material = material;
            end
            
            % add element manually
            obj.Elements(end+1).Element = material;

            obj.Elements(end).IndicesInt = [];
            obj.Elements(end).IndicesExt = material.Space.gnum';
            obj.Elements(end).Indices = [obj.Elements(end).IndicesExt, obj.Elements(end).IndicesInt];
        end

        function patchNumbers = getPatchNumbers(obj, regionNumbers)
            patchNumbers = [];
            for reg = regionNumbers
                patchNumbers = [patchNumbers, obj.Regions(reg).Patches]; %#ok
            end
        end

        function bndNumbers = getBoundaryPatchNumbers(obj, boundaryNumbers)
            bndNumbers = [];
            for bnd = boundaryNumbers
                bndNumbers = [bndNumbers, obj.Boundaries(bnd).PatchBoundary]; %#ok
            end
        end

        function dofs = getBoundaryDOFs(obj, boundaryNumbers)
             dofs = [];
             for bnd = boundaryNumbers
                 localDOFS = obj.Spaces.boundary.gnum{bnd};
                 globalDOFS = obj.Spaces.boundary.dofs(localDOFS);
                 dofs = union(dofs, globalDOFS);
             end
         end

        % Sort DOFs such that the boundary discretization is sorted with the geometric sorting
        function dofs = getBoundaryDOFsSorted(obj, boundaryNumbers)
            dofs = [];
            for iBnd = 1:numel(boundaryNumbers)-1
                % first line element
                bndNr1 = boundaryNumbers(iBnd);
                patch1 = obj.PatchBoundaries(bndNr1).patches;
                side1 = obj.PatchBoundaries(bndNr1).faces;
                nrb1 = nrbextract(obj.Geometry(patch1).nurbs, side1);
                localDOFs1 = obj.Spaces.boundary.gnum{bndNr1};
                globalDOFs1 = obj.Spaces.boundary.dofs(localDOFs1);
                % second line element
                bndNr2 = boundaryNumbers(iBnd+1);
                patch2 = obj.PatchBoundaries(bndNr2).patches;
                side2 = obj.PatchBoundaries(bndNr2).faces;
                nrb2 = nrbextract(obj.Geometry(patch2).nurbs, side2);
                localDOFs2 = obj.Spaces.boundary.gnum{bndNr2};
                globalDOFs2 = obj.Spaces.boundary.dofs(localDOFs2);
                % check if mapping is reversed
                tol = 1e-10;
                if abs(vecnorm(nrbeval(nrb1, 1) - nrbeval(nrb2, 0))) < tol || abs(vecnorm(nrbeval(nrb1, 1) - nrbeval(nrb2, 1))) < tol
                    dofs = union(dofs, globalDOFs1, 'stable');
                elseif abs(vecnorm(nrbeval(nrb1, 0) - nrbeval(nrb2, 0))) < tol || abs(vecnorm(nrbeval(nrb1, 0) - nrbeval(nrb2, 1))) < tol
                    dofs = union(dofs, flip(globalDOFs1), 'stable');
                else
                    error("Boundaries are not matching or numbers not given sorted!")
                end
            end
            localDOFs1 = obj.Spaces.boundary.gnum{boundaryNumbers(end)};
            globalDOFs1 = obj.Spaces.boundary.dofs(localDOFs1);
            % Check last line segment
            if abs(vecnorm(nrbeval(nrb2, 0) - nrbeval(nrb1, 0))) < tol || abs(vecnorm(nrbeval(nrb2, 0) - nrbeval(nrb1, 1))) < tol
                dofs = union(dofs, globalDOFs1, 'stable');
            elseif abs(vecnorm(nrbeval(nrb2, 1) - nrbeval(nrb1, 0))) < tol || abs(vecnorm(nrbeval(nrb2, 1) - nrbeval(nrb1, 1))) < tol
                dofs = union(dofs, flip(globalDOFs1), 'stable');
            else
                error("Boundaries are not matching or numbers not given sorted!")
            end
        end

        function [leftDofs, rightDofs] = getBoundaryDOFsperiodic(obj, boundaryNumbersLeft, boundaryNumbersRight)
            assert(numel(boundaryNumbersLeft) == numel(boundaryNumbersRight), "Number of left and right sides must be the same!");
            leftDofs = [];
            rightDofs = [];
            for iBnd = 1:numel(boundaryNumbersRight)
                % extract left boundary
                bndNrLeft = boundaryNumbersLeft(iBnd);
                patchLeft = obj.PatchBoundaries(bndNrLeft).patches;
                sideLeft = obj.PatchBoundaries(bndNrLeft).faces;
                localDOFsLeft = obj.Spaces.boundary.gnum{bndNrLeft};
                globalDOFsLeft = obj.Spaces.boundary.dofs(localDOFsLeft);
                nrbLeft = nrbextract(obj.Geometry(patchLeft).nurbs, sideLeft);
                % extract right boundary
                bndNrRight = boundaryNumbersRight(iBnd);
                patchRight = obj.PatchBoundaries(bndNrRight).patches;
                sideRight = obj.PatchBoundaries(bndNrRight).faces;
                localDOFsRight = obj.Spaces.boundary.gnum{bndNrRight};
                globalDOFsRight = obj.Spaces.boundary.dofs(localDOFsRight);
                nrbRight = nrbextract(obj.Geometry(patchRight).nurbs, sideRight);
                % check if mapping is reversed
                tol = 1e-10;
                if abs(vecnorm(nrbeval(nrbRight, 0)) - vecnorm(nrbeval(nrbLeft, 0))) < tol
                    leftDofs = union(leftDofs, globalDOFsLeft, 'stable');
                    rightDofs = union(rightDofs, globalDOFsRight, 'stable');
                elseif abs(vecnorm(nrbeval(nrbRight, 0)) - vecnorm(nrbeval(nrbLeft, 1))) < tol
                    leftDofs = union(leftDofs, globalDOFsLeft, 'stable');
                    rightDofs = union(rightDofs, flip(globalDOFsRight), 'stable');
                else
                    error("Boundaries are not matching!")
                end
            end
        end

        %%%%%%%%%%%%%%%% misc
        function calculateArea(obj)
            % patches area and inertia
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).Area = sum(obj.MeshesEval{iPatch}.jacdet .* obj.MeshesEval{iPatch}.quad_weights, 'all');
                % obj.Geometry(iPatch).Inertia = sum(obj.MeshesEval{iPatch}.jacdet .* obj.MeshesEval{iPatch}.quad_weights .* reshape((obj.MeshesEval{iPatch}.geo_map(1, :, :).^2 + obj.MeshesEval{iPatch}.geo_map(2, :, :).^2) , obj.MeshesEval{iPatch}.nqn, obj.MeshesEval{iPatch}.nel), 'all');
            end
            % interfaces lengths
            for iInt = 1:numel(obj.PatchInterfaces)
                patch1 = obj.PatchInterfaces(iInt).patch1;
                side1 = obj.PatchInterfaces(iInt).side1;
                patch2 = obj.PatchInterfaces(iInt).patch2;
                side2 = obj.PatchInterfaces(iInt).side2;
                nrb = obj.Geometry(patch1).boundary(side1).nurbs;
                rule      = msh_gauss_nodes (nrb.order);
                [qn, qw]  = msh_set_quad_nodes (nrb.knots, rule);
                mesh = msh_cartesian(nrb.knots, qn, qw, geo_load(nrb));
                obj.PatchInterfaces(iInt).length = int_Omega(msh_precompute(mesh));
                obj.Geometry(patch1).Lengths{side1} = obj.PatchInterfaces(iInt).length;
                obj.Geometry(patch2).Lengths{side2} = obj.PatchInterfaces(iInt).length;
            end
            % boundary lengths
            for iBnd = 1:numel(obj.PatchBoundaries)
                patch = obj.PatchBoundaries(iBnd).patches;
                side = obj.PatchBoundaries(iBnd).faces;
                nrb = obj.Geometry(patch).boundary(side).nurbs;
                rule      = msh_gauss_nodes (nrb.order);
                [qn, qw]  = msh_set_quad_nodes (nrb.knots, rule);
                mesh = msh_cartesian(nrb.knots, qn, qw, geo_load(nrb));
                obj.PatchBoundaries(iBnd).length = int_Omega(msh_precompute(mesh));
                obj.Geometry(patch).Lengths{side} = obj.PatchBoundaries(iBnd).length;
            end
        end

        function checkPlotValues(obj)
            if obj.PlotValuesNeedUpdate
                obj.calculatePlotValues();
                obj.calculatePlotValuesVector();
            end
            obj.PlotValuesNeedUpdate = false;
        end

        function calculatePlotPoints(obj, nsubs, patches)
            arguments
                obj 
                nsubs = obj.PlotResolution;
                patches = 1:obj.NumberPatches;
            end
            obj.PlotResolution = nsubs;
            obj.PlotValuesNeedUpdate = true;

            for iPatch = patches
                pointsX = unique(kntrefine(obj.Geometry(iPatch).nurbs.knots{1}, obj.PlotResolution(1), obj.Geometry(iPatch).nurbs.order(1), obj.Geometry(iPatch).nurbs.order(1)-1));
                pointsY = unique(kntrefine(obj.Geometry(iPatch).nurbs.knots{2}, obj.PlotResolution(2), obj.Geometry(iPatch).nurbs.order(2), obj.Geometry(iPatch).nurbs.order(2)-1));
                plotPointsParam = {pointsX, pointsY};
                obj.Geometry(iPatch).PlotPointsParam = plotPointsParam;
                obj.Geometry(iPatch).PlotDimensions = [numel(pointsX), numel(pointsY)];
                pnts = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, plotPointsParam, obj.PlotAngle);
                obj.Geometry(iPatch).PlotPoints = pnts(:, 1);
                % update interfaces
                for iInt = 1:numel(obj.PatchInterfaces)
                    if obj.PatchInterfaces(iInt).patch1 == iPatch
                        obj.PatchInterfaces(iInt).PlotPoints = pnts(:, obj.PatchInterfaces(iInt).side1+1);
                        obj.PatchInterfaces(iInt).PlotPointsX = [obj.PatchInterfaces(iInt).PlotPoints{1}(1:end-1)', obj.PatchInterfaces(iInt).PlotPoints{1}(2:end)'];
                        obj.PatchInterfaces(iInt).PlotPointsY = [obj.PatchInterfaces(iInt).PlotPoints{2}(1:end-1)', obj.PatchInterfaces(iInt).PlotPoints{2}(2:end)'];
                        obj.PatchInterfaces(iInt).PlotDimensions = size(obj.PatchInterfaces(iInt).PlotPointsX);
                    end
                end
                % update boundaries
                for iBnd = 1:numel(obj.PatchBoundaries)
                    if obj.PatchBoundaries(iBnd).patches == iPatch
                        obj.PatchBoundaries(iBnd).PlotPoints = pnts(:, obj.PatchBoundaries(iBnd).faces+1);
                        obj.PatchBoundaries(iBnd).PlotPointsX = [obj.PatchBoundaries(iBnd).PlotPoints{1}(1:end-1)', obj.PatchBoundaries(iBnd).PlotPoints{1}(2:end)'];
                        obj.PatchBoundaries(iBnd).PlotPointsY = [obj.PatchBoundaries(iBnd).PlotPoints{2}(1:end-1)', obj.PatchBoundaries(iBnd).PlotPoints{2}(2:end)'];
                        obj.PatchBoundaries(iBnd).PlotDimensions = size(obj.PatchBoundaries(iBnd).PlotPointsX);
                    end
                end
            end
        end

        function updatePlotPoints(obj, rot)
            for iPatch = 1:obj.NumberPatches
                pnts = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsParam, rot);
                obj.Geometry(iPatch).PlotPoints = pnts(:, 1);
                obj.Geometry(iPatch).PlotPointsVector = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsVectorParam, rot);

                % update interfaces
                for iInt = 1:numel(obj.PatchInterfaces)
                    if obj.PatchInterfaces(iInt).patch1 == iPatch
                        obj.PatchInterfaces(iInt).PlotPoints = pnts(:, obj.PatchInterfaces(iInt).side1+1);
                        obj.PatchInterfaces(iInt).PlotPointsX = [obj.PatchInterfaces(iInt).PlotPoints{1}(1:end-1)', obj.PatchInterfaces(iInt).PlotPoints{1}(2:end)'];
                        obj.PatchInterfaces(iInt).PlotPointsY = [obj.PatchInterfaces(iInt).PlotPoints{2}(1:end-1)', obj.PatchInterfaces(iInt).PlotPoints{2}(2:end)'];
                    end
                end
                % update boundaries
                for iBnd = 1:numel(obj.PatchBoundaries)
                    if obj.PatchBoundaries(iBnd).patches == iPatch
                        obj.PatchBoundaries(iBnd).PlotPoints = pnts(:, obj.PatchBoundaries(iBnd).faces+1);
                        obj.PatchBoundaries(iBnd).PlotPointsX = [obj.PatchBoundaries(iBnd).PlotPoints{1}(1:end-1)', obj.PatchBoundaries(iBnd).PlotPoints{1}(2:end)'];
                        obj.PatchBoundaries(iBnd).PlotPointsY = [obj.PatchBoundaries(iBnd).PlotPoints{2}(1:end-1)', obj.PatchBoundaries(iBnd).PlotPoints{2}(2:end)'];
                    end
                end
            end
        end

        function calculatePlotPointsVector(obj, nVectors)
            arguments
                obj 
                nVectors = 300;
            end
            modelArea = sum(vertcat(obj.Geometry.Area));
            for iPatch = 1:obj.NumberPatches
                L1 = obj.Geometry(iPatch).Lengths{1} + obj.Geometry(iPatch).Lengths{2};
                L2 = obj.Geometry(iPatch).Lengths{3} + obj.Geometry(iPatch).Lengths{4};
                num_vectors_i = ceil(obj.Geometry(iPatch).Area/modelArea*nVectors);
                % Number vectors in x an y according to lengths of sides
                q1 = ceil((num_vectors_i* L2/L1)^0.5);
                q2 = ceil((num_vectors_i* L1/L2)^0.5);
                % Points with also half distance to boundary
                paramPlotPoints = {linspace(0, 1-1/q1, q1)+0.5/q1, linspace(0, 1-1/q2, q2)+0.5/q2};
                % Points also directly on boundary
                %paramPlotPoints = {linspace(0, 1, q1), linspace(0, 1, q2)};
                obj.Geometry(iPatch).PlotPointsVectorParam = paramPlotPoints;
                obj.Geometry(iPatch).PlotPointsVector = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsVectorParam, obj.PlotAngle);
                obj.Geometry(iPatch).PlotDimensionsVector = size(obj.Geometry(iPatch).PlotPointsVector{1,1});
            end
        end

        function calculatePlotValuesVector(obj)
            progressbar([obj.Name ': Calculating function space at vector plot points: '])
            for iPatch = 1:obj.NumberPatches
                [obj.Geometry(iPatch).PlotMeshVector, obj.Geometry(iPatch).PlotSpaceVector] = sp_eval_pnts(obj.Spaces.sp_patch{iPatch}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsVectorParam); 
                progressbar(iPatch/obj.NumberPatches);
            end
            progressbar('finished');
        end

        function setEvaluationPoints(obj, xPhysical)
            arguments
                obj IGA_Region;
                xPhysical (2,:) double
            end
            [xParam, patches] = physical2param(vertcat(obj.Geometry.nurbs), xPhysical);
            obj.setEvaluationPointsParam(xParam, patches);
        end

        function setEvaluationPointsParam(obj, xParam, patchNumbers)
            arguments
                obj IGA_Region;
                xParam (2,:) double
                patchNumbers (1,:) double
            end
            obj.EvaluationPointsParam = xParam;
            obj.EvaluationPatches = patchNumbers;
            obj.EvaluationShapeValues = shp_eval_mp(obj.Spaces, obj.Geometry, xParam, patchNumbers, {'value', 'gradient'});
        end
        
        function [G, H] = getBCmatrices(obj)
            [G, H] = getBCmatrices@EQS_Element(obj);
        end

        function [dofsR, dofsK] = getBCindices(obj)
            [dofsR, dofsK] = getBCindices@EQS_Element(obj);
        end

        function [b, db] = getBCvalues(obj, t)
            [b, db] = getBCvalues@EQS_Element(obj, t);
        end

        function W = OmegaMatrix(obj)
            W = obj.omega * speye(obj.NumberDOF);
        end

        %%%%%%%%%%%%%%%% plots
        function plotGeometry(obj, t, args)
            arguments
                obj
                t = 0;
                args.FocusFigure (1,1) logical = true;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
            end
            
            obj.plotRegions(t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", args.StoreName, "FocusFigure", args.FocusFigure);
            obj.plotBoundary(t, "FigureName", args.FigureName, "ValueName", args.ValueName, "FocusFigure", args.FocusFigure);
            daspect([1,1,1]);
            % axis equal
        end

        function plotRegions(obj, t, args)
            arguments
                obj
                t = 0;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkRotationUpdate(t);

            obj.PlotTool.initializeColor(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensions), args.FocusFigure);
            
            for iRegion = 1:obj.NumberRegions
                patches = obj.Regions(iRegion).Patches;
                if isempty(obj.Regions(iRegion).Material)
                    col = "white";
                    alph = 1;
                    warning("IGA_Region:MaterialNotDefined: Material for region %i is not yet defined!", iRegion)
                else
                    col = obj.Regions(iRegion).Material.getPlotColor();
                    alph = obj.Regions(iRegion).Material.getPlotAlpha();
                end
                for patch = patches
                    obj.PlotTool.setColor(args.FigureName, args.ValueName, args.StoreName, patch, "XData", obj.Geometry(patch).PlotPoints{1}, ...
                        "YData", obj.Geometry(patch).PlotPoints{2}, "Color", col, "Alpha", alph);
                end
            end 
            drawOrder = horzcat(obj.Regions.Patches);
            obj.PlotTool.updateColor(args.FigureName, args.ValueName, args.StoreName, drawOrder);
        end

        function plotBoundary(obj, t, args)
            arguments
                obj 
                t = 0;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.LineWidthMaterials = 1.5;
                args.ColorMaterials = "black";
                args.LineWidthInterfaces = 0.5;
                args.ColorInterfaces = "none";
                args.LineWidthBoundaries = 1.5;
                args.ColorBoundaries = "black";
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkRotationUpdate(t);

            if ~isempty(obj.PatchInterfaces)
                % plot patch interfaces with no material jump
                InterfacesNoMat = obj.PatchInterfaces(vertcat(obj.PatchInterfaces.IsMaterialInterface) == 0);
                if ~isempty(InterfacesNoMat)
                    obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "InterfacesNoMat", vertcat(InterfacesNoMat.PlotDimensions), args.FocusFigure);
                    for iInt = 1:numel(InterfacesNoMat)
                        obj.PlotTool.setLines(args.FigureName, args.ValueName, "InterfacesNoMat", iInt, ...
                            "XData", InterfacesNoMat(iInt).PlotPointsX, "YData", InterfacesNoMat(iInt).PlotPointsY);
                    end
                    obj.PlotTool.updateLines(args.FigureName, args.ValueName, "InterfacesNoMat", "LineWidth", args.LineWidthInterfaces, "EdgeColor", args.ColorInterfaces);
                end
    
                % plot patch interfaces with material jump
                InterfacesMat = obj.PatchInterfaces(vertcat(obj.PatchInterfaces.IsMaterialInterface) == 1);
                if ~isempty(InterfacesMat)
                    obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "InterfacesMat", vertcat(InterfacesMat.PlotDimensions), args.FocusFigure);
                    for iInt = 1:numel(InterfacesMat)
                        obj.PlotTool.setLines(args.FigureName, args.ValueName, "InterfacesMat", iInt, ...
                            "XData", InterfacesMat(iInt).PlotPointsX, "YData", InterfacesMat(iInt).PlotPointsY);
                    end
                    obj.PlotTool.updateLines(args.FigureName, args.ValueName, "InterfacesMat", "LineWidth", args.LineWidthMaterials, "EdgeColor", args.ColorMaterials);
                end
            end

            % plot boundaries
            obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "Boundaries", vertcat(obj.PatchBoundaries.PlotDimensions), args.FocusFigure);
            for iBnd = 1:numel(obj.PatchBoundaries)
                obj.PlotTool.setLines(args.FigureName, args.ValueName, "Boundaries", iBnd, ...
                    "XData", obj.PatchBoundaries(iBnd).PlotPointsX, "YData", obj.PatchBoundaries(iBnd).PlotPointsY);
            end
            obj.PlotTool.updateLines(args.FigureName, args.ValueName, "Boundaries", "LineWidth", args.LineWidthBoundaries, "EdgeColor", args.ColorBoundaries);
        end

        function plotKnotLines(obj)
            for iPatch = 1:obj.NumberPatches
                msh = obj.Meshes.msh_patch{iPatch};
                for knt = unique(obj.Spaces.sp_patch{iPatch}.knots{2})
                    knot_vals =  linspace(0,1,20);
                    xy = msh.map({knot_vals, knt});
                    line(xy(1,:), xy(2,:), 'Color', [1,1,1]*0.5);
                end
                % Y lines
                for knt = unique(obj.Spaces.sp_patch{iPatch}.knots{1})
                    knot_vals =  linspace(0,1,20);
                    xy = msh.map({knt, knot_vals});
                    line(xy(1,:), xy(2,:), 'Color', [1,1,1]*0.5);
                end
            end
        end
        
        function exportIGES1D(obj, args)
            arguments
                obj 
                args.Filename = 'nrbexport.igs';
                args.Scaling = 1;
                args.IncludePatchInterfaces (1,1) logical = false;
            end

            % Add all NURBS boundaries to export
            for iBnd = 1:numel(obj.PatchBoundaries)
                patch = obj.PatchBoundaries(iBnd).patches;
                side = obj.PatchBoundaries(iBnd).faces;
                nrb_export(iBnd) = nrbextract(obj.Geometry(patch).nurbs, side);
            end
            % Add NURBS interfaces to export
            for iInt = 1:numel(obj.PatchInterfaces)
                if obj.PatchInterfaces(iInt).IsMaterialInterface || args.IncludePatchInterfaces
                    patch = obj.PatchInterfaces(iInt).patch1;
                    side = obj.PatchInterfaces(iInt).side1;
                    nrb_export(end+1) = nrbextract(obj.Geometry(patch).nurbs, side);
                end
            end
            % Scale
            for iNrb = 1:numel(nrb_export)
                nrb_export(iNrb).coefs(1:3,:,:) = nrb_export(iNrb).coefs(1:3,:,:)*args.Scaling;
            end
            % Export to IGES
            nrb2iges(nrb_export, args.Filename)
        end

        function exportIGES2D(obj, args)
            arguments
                obj 
                args.Filename = 'nrbexport.igs';
                args.Scaling = 1;
            end
            nrb_export = vertcat(obj.Geometry.nurbs);
            for iGeo = 1:numel(obj.Geometry)
                nrb_export(iGeo).coefs(1:3,:,:) = nrb_export(iGeo).coefs(1:3,:,:)*args.Scaling;
            end
            nrb2iges(nrb_export, args.Filename)
        end

        function plotBoundaryConditions(obj)
            lw = 2;
            for iBC = 1:numel(obj.BoundaryConditions)
                checkBC = obj.BoundaryConditions(iBC).BoundaryCondition;
                if isprop(checkBC, "BoundaryPatchNumbers")
                    col = checkBC.getPlotColor();
                    for bndNr =  checkBC.BoundaryPatchNumbers
                        plot(obj.PatchBoundaries(bndNr).PlotPoints{1}, obj.PatchBoundaries(bndNr).PlotPoints{2}, "Color", col, "LineWidth", lw);
                    end
                end
            end
        end

        function plotPatchNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14;
                args.Color = "black";
            end
            xmid = zeros(2, obj.NumberPatches);
            for iPatch = 1:obj.NumberPatches
                xmid(:,iPatch) = obj.Meshes.msh_patch{iPatch}.map([0.5, 0.5]);
            end
            text(xmid(1,:), xmid(2,:), string(1:obj.NumberPatches), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize, "Color", args.Color);
        end

        function plotRegionNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14;
                args.Color = "black";
                args.RegionNumbers (1,:) double = 1:obj.NumberRegions;
                args.SnapNearest (1,1) logical = true;  % false center (can be outside region) or snap to nearest patch center
            end
            obj.checkRotationUpdate();
            xmidRegions = zeros(2, numel(args.RegionNumbers));
            for iRegion = 1:numel(args.RegionNumbers)
                xmidRegions(:, iRegion) = obj.getRegionCenter(args.RegionNumbers(iRegion), args.SnapNearest);
            end
            text(xmidRegions(1,:), xmidRegions(2,:), string(args.RegionNumbers), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize, "Color", args.Color);
        end

        function xCenter = getRegionCenter(obj, regionNumbers, snapNearest)
            arguments
                obj 
                regionNumbers 
                snapNearest (1,1) logical = false;  % true center (can be outside region) or snap to nearest patch center
            end
            patches = obj.getPatchNumbers(regionNumbers);
            xMidPatches = zeros(3, numel(patches));
            for iPatch = 1:numel(patches)
                xMidPatches(:, iPatch) = nrbeval(obj.Geometry(patches(iPatch)).nurbs, {0.5, 0.5});
            end
            AreaMoment = [0;0];
            Area = 0;
            for iPatch = 1:numel(patches)
                AreaMoment = AreaMoment + obj.Geometry(patches(iPatch)).Area*xMidPatches(1:2, iPatch);
                Area = Area + obj.Geometry(patches(iPatch)).Area;
            end
            xMean = AreaMoment/Area;
            if snapNearest
                [~, ind] = min(vecnorm(xMean - xMidPatches(1:2,:)));
                xCenter = xMidPatches(1:2, ind);
            else
                xCenter = xMean;
            end
            xCenter = [cos(obj.RotationAngle), -sin(obj.RotationAngle); sin(obj.RotationAngle), cos(obj.RotationAngle)]*xCenter;
        end

        function plotBoundaryNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14;
                args.Color = "black";
            end
            xmidPatchBoundaries = zeros(3, obj.NumberBoundaries);
            for iBndPatch = 1:numel(obj.PatchBoundaries)
                xmidPatchBoundaries(:,iBndPatch) = nrbeval(obj.Geometry(obj.PatchBoundaries(iBndPatch).patches).boundary(obj.PatchBoundaries(iBndPatch).faces).nurbs, 0.5);
            end
            xmidBoundaries = zeros(3, obj.NumberBoundaries);
            for iBnd = 1:obj.NumberBoundaries
                Lengths = 0;
                xcent = [0;0;0];
                for iBndPatch = 1:numel(obj.Boundaries(iBnd).PatchBoundary)
                    patchBoundary = obj.Boundaries(iBnd).PatchBoundary(iBndPatch);
                    lBnd = obj.PatchBoundaries(patchBoundary).length;
                    xcent = xcent + xmidPatchBoundaries(:, patchBoundary)*lBnd;
                    Lengths = Lengths + lBnd;
                end
                xmid = xcent / Lengths;
                xmidtest = xmidPatchBoundaries(:, obj.Boundaries(iBnd).PatchBoundary);
                [~, ind] = min(vecnorm(xmid - xmidtest));
                xmidBoundaries(:, iBnd) = xmidtest(:, ind);
            end
            text(xmidBoundaries(1,:), xmidBoundaries(2,:), string(1:obj.NumberBoundaries), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize, "Color", args.Color);
        end

        function plotPatchBoundaryNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14; 
            end
            for iBnd = 1:numel(obj.PatchBoundaries)
                x = nrbeval(obj.Geometry(obj.PatchBoundaries(iBnd).patches).boundary(obj.PatchBoundaries(iBnd).faces).nurbs, 0.5);
                text(x(1), x(2), num2str(iBnd), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize);
            end
        end

        function plotPatchInterfaceNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14; 
            end
            for iInt = 1:numel(obj.PatchInterfaces)
                x = nrbeval(obj.Geometry(obj.PatchInterfaces(iInt).patch1).boundary(obj.PatchInterfaces(iInt).side1).nurbs, 0.5);
                text(x(1), x(2), num2str(iInt), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize);
            end
        end

        function plotPatchBoundaryInterfaceSides(obj, bndNr)
            patch = obj.PatchBoundaries(bndNr).patches;
            side = obj.PatchBoundaries(bndNr).faces;
            x1 = nrbeval(obj.Geometry(patch).boundary(side).nurbs, 0);
            x2 = nrbeval(obj.Geometry(patch).boundary(side).nurbs, 1);
            text(x1(1), x1(2), "1", "HorizontalAlignment", "center", "VerticalAlignment", "middle", "Color", "red", "FontSize", 14);
            text(x2(1), x2(2), "2", "HorizontalAlignment", "center", "VerticalAlignment", "middle", "Color", "red", "FontSize", 14);
        end

        function moveInterfaceToBoundary(obj, interfaceNr)
            for iInt = 1:numel(interfaceNr)
                obj.PatchBoundaries(end+1).nsides = 1;
                obj.PatchBoundaries(end).patches = obj.PatchInterfaces(interfaceNr(iInt)).patch1;
                obj.PatchBoundaries(end).faces = obj.PatchInterfaces(interfaceNr(iInt)).side1;

                obj.PatchBoundaries(end+1).nsides = 1;
                obj.PatchBoundaries(end).patches = obj.PatchInterfaces(interfaceNr(iInt)).patch2;
                obj.PatchBoundaries(end).faces = obj.PatchInterfaces(interfaceNr(iInt)).side2;
            end
            obj.PatchInterfaces(interfaceNr) = [];
            warning("BoundaryInterfaces might not be correct, you need to add them manually!")
        end

        function boundaryNumbers = getBoundaryNumbers(obj, patchNumbers, patchSides)
            % if only the struct with patches and sides is given: extract them
            if nargin == 2
                patchSides = patchNumbers.sides;
                patchNumbers = patchNumbers.patches;
            end

            boundaryNumbers = [];
            for iPatch = 1:numel(patchNumbers)
                for iBnd = 1:numel(obj.PatchBoundaries)
                    if patchNumbers(iPatch) == obj.PatchBoundaries(iBnd).patches && patchSides(iPatch) == obj.PatchBoundaries(iBnd).faces
                        boundaryNumbers = [boundaryNumbers, iBnd];
                        break;
                    end
                    if iBnd == numel(obj.PatchBoundaries)
                        warning('IGA_Region:boundaryNumbers', 'Some boundary numbers were not found! \t Patch: %i \t Side: %i \n', patchNumbers(iPatch) , patchSides(iPatch));
                    end
                end
            end
        end

        function plotControlPoints(obj, args)
            arguments
                obj 
                args.DrawText (1,1) logical = true;
                args.FontSize = 14;
                args.Indices = 1:size(obj.ControlPoints, 1);
                args.MarkerSize = 15;
                args.Color = TUDa_getColor("9b");
            end

            scatter(obj.ControlPoints(args.Indices, 1), obj.ControlPoints(args.Indices, 2), "filled", 'MarkerFaceColor', args.Color, 'SizeData', args.MarkerSize, Marker='o')
            if args.DrawText
                text(obj.ControlPoints(args.Indices, 1), obj.ControlPoints(args.Indices, 2), string(args.Indices), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize);
            end
        end

        function plotQuadraturePoints(obj, patches, markerSize, col)
            if ~exist("patches", "var") || isempty(patches)
                patches = 1:obj.NumberPatches;
            end
            if ~exist("markerSize", "var") || isempty(markerSize)
                markerSize = 15;
            end
            if ~exist("col", "var")
                col = TUDa_getColor("9b");
            end
            points = [];
            for patch = patches
                points = [points, reshape(obj.MeshesEval{patch}.geo_map, 2, [])];

            end
            scatter(points(1,:), points(2,:), "filled", 'MarkerFaceColor', col, 'SizeData', markerSize, Marker='o');
        end

        function plotEvaluationPoints(obj, options)
            arguments
                obj = IGA_Region;
                options.PlotNumbers (1,1) logical = false;
                options.FontSize = 16;
            end
            xPhysical = zeros(size(obj.EvaluationPointsParam));
            for iPoint = 1:size(obj.EvaluationPointsParam, 2)
                xPhysical(:,iPoint) = obj.Meshes.msh_patch{obj.EvaluationPatches(iPoint)}.map(obj.EvaluationPointsParam(:, iPoint));
            end
            scatter(xPhysical(1,:), xPhysical(2,:), "black", "filled");
            if options.PlotNumbers
                text(xPhysical(1,:), xPhysical(2,:), " " +string(1:size(xPhysical,2)), "HorizontalAlignment", "left", "VerticalAlignment", "middle", "FontSize", options.FontSize);
            end
        end

        % Calculate the Shape Values and Gradients at Plot Points in advance
        function calculatePlotValues(obj)
            progressbar([obj.Name ': Calculating function space at plot points: '])
            for iPatch = 1:obj.NumberPatches
                [obj.Geometry(iPatch).PlotMeshScalar, obj.Geometry(iPatch).PlotSpaceScalar] = sp_eval_pnts(obj.Spaces.sp_patch{iPatch}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsParam); 
                progressbar(iPatch/obj.NumberPatches);
            end
            progressbar('finished');
        end

        function pnts = getGeoPlotPoints(~, nurbs, pts, rot)
            if nargin <= 3
                rot = 0;
            end
            xpts = pts{1};
            ypts = pts{2};

            F_geo = nrbeval(nurbs, {xpts, ypts});
            X = squeeze(F_geo(1, :, :))*cos(rot) - squeeze(F_geo(2, :, :))*sin(rot);
            Y = squeeze(F_geo(2, :, :))*cos(rot) + squeeze(F_geo(1, :, :))*sin(rot);
            pnts{1, 1} = X;
            pnts{2, 1} = Y;

            bnd = nrbeval(nrbextract(nurbs, 1), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 2} = X;
            pnts{2, 2} = Y;

            bnd = nrbeval(nrbextract(nurbs, 2), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 3} = X;
            pnts{2, 3} = Y;

            bnd = nrbeval(nrbextract(nurbs, 3), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 4} = X;
            pnts{2, 4} = Y;

            bnd = nrbeval(nrbextract(nurbs, 4), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 5} = X;
            pnts{2, 5} = Y;
        end

        function exportXML(obj, name)
            disp("THIS IS WORK IN PROGRESS: connectivity and nurbs weights not yet included")
            fileID = fopen(name, 'w');
            fprintf(fileID, '<xml>\n');
            arr2text = @(text) regexprep(num2str(text),'\s+',' ');
            for iPatch = 1:obj.NumberPatches
                fprintf(fileID, '  <!-- Patch %i -->\n', iPatch);
                type = 'TensorBSpline2';
                id = iPatch;
                basis2d = 'TensorBSplineBasis2';
                basis1d = 'BSplineBasis';
                fprintf(fileID, '  <Geometry type="%s" id="%i">\n', type, id);
                fprintf(fileID, '    <Basis type="%s">\n', basis2d);
                fprintf(fileID, '      <Basis type="%s" index="0">\n', basis1d);
                fprintf(fileID, '        <KnotVector degree="%i">%s </KnotVector>\n', obj.Geometry(iPatch).nurbs.order(1)-1, arr2text(obj.Geometry(iPatch).nurbs.knots{1}));
                fprintf(fileID, '      </Basis>\n');
                fprintf(fileID, '      <Basis type="%s" index="1">\n', basis1d);
                fprintf(fileID, '        <KnotVector degree="%i">%s </KnotVector>\n', obj.Geometry(iPatch).nurbs.order(2)-1, arr2text(obj.Geometry(iPatch).nurbs.knots{2}));
                fprintf(fileID, '      </Basis>\n');
                fprintf(fileID, '    </Basis>\n');
                fprintf(fileID, '    <coefs geoDim="%i">\n', obj.Geometry(iPatch).rdim);
                coefs = reshape(obj.Geometry(iPatch).nurbs.coefs(1:obj.Geometry(iPatch).rdim,:,:), obj.Geometry(iPatch).rdim, [])';
                for iCtrlpnt = 1:size(coefs, 1)
                    fprintf(fileID, '      %s\n', arr2text(coefs(iCtrlpnt, :)));
                end
                fprintf(fileID, '    </coefs>\n');
                fprintf(fileID, '  </Geometry>\n');
                fprintf(fileID, '\n');
            end
            % interfaces
            fprintf(fileID, '  <!-- Multipatch -->\n');
            parDim = 2;
            id = 0;
            fprintf(fileID, '  <MultiPatch parDim="%i" id="%i">\n', parDim, id);
            fprintf(fileID, '    <patches type="id_range">%i %i</patches>\n', 1, obj.NumberPatches);
            fprintf(fileID, '    <interfaces>\n');
            for iInt = 1:numel(obj.PatchInterfaces)
                interfaces = [obj.PatchInterfaces(iInt).patch1, obj.PatchInterfaces(iInt).side1, obj.PatchInterfaces(iInt).patch2, obj.PatchInterfaces(iInt).side2, [0,1,1,1]];
                fprintf(fileID, '      %s \n', arr2text(interfaces)); %1 2 0 1 0 1 1 1
            end

            fprintf(fileID, '    </interfaces>\n');

            fprintf(fileID, '    <boundary>\n');
            for iBnd = 1:numel(obj.PatchBoundaries)
                boundary = [obj.PatchBoundaries(iBnd).patches, obj.PatchBoundaries(iBnd).faces];
                fprintf(fileID, '      %s \n', arr2text(boundary));
            end
            fprintf(fileID, '    </boundary>\n');
            fprintf(fileID, '  </MultiPatch>\n');

            fprintf(fileID, '</xml>');
            fclose(fileID);
        end
    
        function checkRotationUpdate(obj, t)
            if exist("t", "var") && isfield(obj.Solution, "RotationAngle") && ~isscalar(obj.Solution.RotationAngle)
                obj.RotationAngle = interp1(obj.Solution.Time, obj.Solution.RotationAngle, t);
            end
            if obj.RotationAngle ~= obj.PlotAngle
                obj.updatePlotPoints(obj.RotationAngle);
                obj.PlotAngle = obj.RotationAngle;
            end
        end
    end
end
