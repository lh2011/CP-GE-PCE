classdef FEM_Region < EQS_System

    properties
        DegreeQuadrature;
        % domain definitions
        Length;                             % axial length
        Geometry;                           % parametrization of the mapping to the physical domain and derivatives
        Boundaries = struct("lines", {}, "dofs", {});                         % information about the patches with an outer boundary
        Triangles; Nodes;
        NodesPlot;              % (rotated) Nodes for plotting

        % computation of spaces and meshes
        Spaces, SpacesEval;                 % space stores the shape functions of discrete spaces
        Meshes, MeshesEval;                 % mesh stores partition of domain and quadrature rule
        NumberRegions;

        Regions = struct('Name', {}, 'SubRegions', {}, 'Material', {});

        % plotting
        UpdateGeometryPlot = true;          % check if plot needs update

        PlotTool;
        RotationAngle = 0;      % current rotation angle of geometry
        PlotAngle = 0;          % angle when geometry was last plotted, if ~=RotationAngle, update plot points
    end

    methods (Abstract)
        generateSpaces(obj)
    end

    methods
        function obj = FEM_Region(length, args)
            arguments
                length (1,1) double = 1;
                args.degreeQuadrature = 1;
            end
            obj.Length = length;
            obj.DegreeQuadrature = args.degreeQuadrature;
            obj.resetMaterials();
            obj.resetBoundaryConditions();
            obj.PlotTool = PLT_TriangularPattern();
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
            
        function generateMeshes(obj, degreeQuadrature, length)
            for iRegion = 1:obj.NumberRegions
                obj.MeshesEval{iRegion} = MESH_2Dtri(obj.Nodes, obj.Geometry(iRegion).Triangles, "DegreeQuadrature", degreeQuadrature, "length", length);
                obj.Geometry(iRegion).Area = int_Omega(obj.MeshesEval{iRegion});
            end
        end

        function importGMSH_m(obj, msh, scaling)
            arguments
                obj;
                msh;
                scaling = 1;
            end
            obj.Nodes = msh.POS(:, 1:2)*scaling;
            obj.Triangles = msh.TRIANGLES;
            % individual regions
            regionNumbers = unique(msh.TRIANGLES(:,4));
            obj.NumberRegions = numel(regionNumbers);
            obj.resetMaterials();
            for iRegion = 1:obj.NumberRegions
                num = regionNumbers(iRegion);
                obj.Geometry(iRegion).Triangles = msh.TRIANGLES(msh.TRIANGLES(:,4) == num, 1:3);
                obj.Regions(iRegion).SubRegions = iRegion;
            end
            obj.NumberDOF = size(obj.Nodes, 1);
            
            lineIDs = unique(msh.LINES(:,3));
            for iLine = 1:numel(lineIDs)
                lineID = lineIDs(iLine);
                obj.Boundaries(iLine).lines = msh.LINES(msh.LINES(:,3)==lineID, 1:2);
                obj.Boundaries(iLine).dofs = unique(obj.Boundaries(iLine).lines);
                obj.Boundaries(iLine).PlotDimensions = size(obj.Boundaries(iLine).lines);
            end

            obj.generateMeshes(obj.DegreeQuadrature, obj.Length);
            obj.generateSpaces();
            obj.NodesPlot = obj.Nodes;
        end

        function importFEMM_ans(obj, ansFile, scaling)
            arguments
                obj 
                ansFile 
                scaling = 1;
            end

            femmdata = read_femmdata(num2str(ansFile));
            obj.Nodes = femmdata.node(:, 1:2)*scaling;
            obj.Triangles = femmdata.elem+1;
            regionNumbers = unique(femmdata.elem(:,4));
            obj.NumberRegions = numel(regionNumbers);
            obj.resetMaterials();
            for iRegion = 1:obj.NumberRegions
                num = regionNumbers(iRegion);
                % +1 as FEMM starts couting from 0
                T = femmdata.elem(femmdata.elem(:,4) == num, 1:3)+1;
                obj.Geometry(iRegion).Triangles = T;
                obj.Regions(iRegion).SubRegions = iRegion;
            end
            % calculate all boundaries
            geometryPoints = femmdata.points*scaling;
            obj.processLinesFEMM(femmdata.segments, geometryPoints);
            obj.processArcsFEMM(femmdata.arcsegments, geometryPoints);

            obj.generateMeshes(obj.DegreeQuadrature, obj.Length);
            obj.generateSpaces();
            obj.NodesPlot = obj.Nodes;
        end

        function processLinesFEMM(obj, segments, points)
            for iSeg = 1:size(segments, 1)
                p1 = points(segments(iSeg,1)+1, 1:2);
                p2 = points(segments(iSeg,2)+1, 1:2);
                v = p2-p1;
                lambda = sum((obj.Nodes-p1).*v, 2)/vecnorm(v)^2;
                normals = obj.Nodes-p1-lambda.*v;
                normals_len = vecnorm(normals,2,2);
                tol = 1e-10;
                onLineInd = find(lambda>-tol & lambda < 1+tol & normals_len<tol);
                lambdaOnLine = lambda(onLineInd);
                [~, indSort] = sort(lambdaOnLine);
                onLineInd =  onLineInd(indSort);

                obj.Boundaries(iSeg).dofs = onLineInd;
                obj.Boundaries(iSeg).lines = [onLineInd(1:end-1), onLineInd(2:end)];
                obj.Boundaries(iSeg).PlotDimensions = size(obj.Boundaries(iSeg).lines);
            end
        end

        function processArcsFEMM(obj, arcs, points)
            nLines = numel(obj.Boundaries);
            for iSeg = 1:size(arcs, 1)
                p1 = points(arcs(iSeg,1)+1, 1:2);
                p2 = points(arcs(iSeg,2)+1, 1:2);
                alpha = arcs(iSeg, 3);
                
                mid = 0.5*(p1+p2);
                secant = p2 - p1;
                nArc = [-secant(2), secant(1)];
                nArc = nArc/vecnorm(nArc);
                dMid = vecnorm(secant);
                R = dMid/(2*sind(alpha/2));
                l = sqrt(R.^2 - (dMid/2).^2);
                center = mid + l*nArc; % Centre of circle point
                % rotate to reference frame
                phi = atan2(-nArc(2), -nArc(1));
                Rmat = [cos(phi), -sin(phi); sin(phi), cos(phi)];
                v1 = (p1 - center)*Rmat;
                v2 = (p2 - center)*Rmat;
                a1 = atan2(v1(2), v1(1));
                a2 = atan2(v2(2), v2(1));

                Rnodes = vecnorm(obj.Nodes - center, 2, 2);
                nodesRotate = (obj.Nodes-center)*Rmat;
                angles = atan2(nodesRotate(:,2), nodesRotate(:,1));

                tol = 1e-10;
                onLineInd = find(angles>a1-tol & angles < a2+tol & abs(R-Rnodes)<tol);
                anglesOnLine = angles(onLineInd);
                [~, indSort] = sort(anglesOnLine);
                onLineInd =  onLineInd(indSort);

                obj.Boundaries(iSeg+nLines).dofs = onLineInd;
                obj.Boundaries(iSeg+nLines).lines = [onLineInd(1:end-1), onLineInd(2:end)];
                obj.Boundaries(iSeg+nLines).PlotDimensions = size(obj.Boundaries(iSeg+nLines).lines);
            end
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
                error("FEM_Region:NotValidRegion: Region number %i is not with %i number of regions", regions, obj.NumberRegions);
            end
            % TBD warning if material already set
            for reg = regions
                obj.Regions(reg).Material = material;
            end
            
            % add element manually
            obj.Elements(end+1).Element = material;

            obj.Elements(end).IndicesInt = [];
            obj.Elements(end).IndicesExt = material.Space.gnum';
            obj.Elements(end).Indices = [obj.Elements(end).IndicesExt, obj.Elements(end).IndicesInt];
        end

        function dofs = getBoundaryDOFs(obj, boundaryNumbers)
            dofs = [];
            for bnd = boundaryNumbers
                dofs = union(obj.Boundaries(bnd).dofs, dofs);
            end
        end

        function [leftDofs, rightDofs] = getBoundaryDOFsperiodic(obj, boundaryNumbersLeft, boundaryNumbersRight)
            assert(numel(boundaryNumbersLeft) == numel(boundaryNumbersRight), "Number of left and right sides must be the same!");
            leftDofs = [];
            rightDofs = [];
            % gather left and right nodes
            for iBnd = 1:numel(boundaryNumbersLeft)
                bndNr = boundaryNumbersLeft(iBnd);
                leftDofs = union(leftDofs, obj.Boundaries(bndNr).dofs);
            end
            for iBnd = 1:numel(boundaryNumbersRight)
                bndNr = boundaryNumbersRight(iBnd);
                rightDofs = union(rightDofs, obj.Boundaries(bndNr).dofs);
            end
            % sort by radial distance
            [~, indLeft] = sort(vecnorm(obj.Nodes(leftDofs, :), 2, 2));
            [~, indRight] = sort(vecnorm(obj.Nodes(rightDofs, :), 2, 2));
            leftDofs = leftDofs(indLeft);
            rightDofs = rightDofs(indRight);
        end

        function plotGeometry(obj, t, args)
            arguments
                obj
                t = 0;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
                args.EdgeColor = "none";
                args.LineWidth = 0.5;
            end
            obj.checkRotationUpdate(t);

            obj.plotRegions(t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName",args.StoreName, "edgeColor", args.EdgeColor, "lineWidth", args.LineWidth);
            obj.plotBoundary(t, args.FigureName, args.ValueName);
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
                args.edgeColor = "none";
                args.lineWidth = 0.5;
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkRotationUpdate(t);

            obj.PlotTool.initializeColor(args.FigureName, args.ValueName, args.StoreName, obj.NumberRegions, args.FocusFigure);

            for iRegion = 1:obj.NumberRegions
                if isempty(obj.Regions(iRegion).Material)
                    col = "white";
                    alph = 1;
                    warning("FEM_Region:MaterialNotDefined: Material for region %i is not yet defined!", iRegion)
                else
                    col = obj.Regions(iRegion).Material.getPlotColor();
                    alph = obj.Regions(iRegion).Material.getPlotAlpha();
                end
                obj.PlotTool.setColor(args.FigureName, args.ValueName, args.StoreName, iRegion, "Faces", obj.Geometry(iRegion).Triangles, ...
                    "Vertices", obj.NodesPlot, "Color", col, "Alpha", alph);
            end
            obj.PlotTool.updateColor(args.FigureName, args.ValueName, args.StoreName, "EdgeColor", args.edgeColor, "LineWidth", args.lineWidth);
        end

        function plotDOFs(obj, dofs, plotText)
            if ~exist("dofs", "var") || isempty(dofs)
                dofs = 1:obj.NumberDOF;
            end
            if ~exist("plotText", "var")
                plotText = false;
            end
            scatter(obj.NodesPlot(dofs, 1), obj.NodesPlot(dofs, 2), "filled", "black");
            if plotText
                text(obj.NodesPlot(dofs, 1), obj.NodesPlot(dofs, 2), string(dofs), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", 16);
            end
        end

        function subRegions = getPatchNumbers(obj, regionNumbers)
            subRegions = [];
            for reg = regionNumbers
                subRegions = [subRegions, obj.Regions(reg).SubRegions]; %#ok
            end
        end

        function plotBoundary(obj, t, FigureName, ValueName, args)
            arguments
                obj 
                t = 0;
                FigureName = "Geometry";
                ValueName = FigureName;
                args.LineWidthBoundaries = 1.5;
                args.ColorBoundaries = "black";
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkRotationUpdate(t);

            % plot boundaries
            obj.PlotTool.initializeLines(FigureName, ValueName, "Boundaries", vertcat(obj.Boundaries.PlotDimensions), args.FocusFigure);
            
            for iBnd = 1:numel(obj.Boundaries)
                xValues = [obj.NodesPlot(obj.Boundaries(iBnd).lines(:,1), 1), obj.NodesPlot(obj.Boundaries(iBnd).lines(:,2), 1)];
                yValues = [obj.NodesPlot(obj.Boundaries(iBnd).lines(:,1), 2), obj.NodesPlot(obj.Boundaries(iBnd).lines(:,2), 2)];
                obj.PlotTool.setLines(FigureName, ValueName, "Boundaries", iBnd, ...
                    "XData", xValues, "YData", yValues);
            end
            obj.PlotTool.updateLines(FigureName, ValueName, "Boundaries", "LineWidth", args.LineWidthBoundaries, "EdgeColor", args.ColorBoundaries);
        end

        function xCenter = getRegionCenter(obj, regionNumber, snapNearest)
            arguments
                obj 
                regionNumber 
                snapNearest (1,1) logical = false;  % true center (can be outside region) or snap to nearest patch center
            end
            
            regionNodeNumbers = unique(obj.Geometry(regionNumber).Triangles);
            regionNodes = obj.Nodes(regionNodeNumbers, :);
            regionMesh = MESH_2D(obj.MeshesEval(regionNumber));
            Area = sum(regionMesh.jacdet.*regionMesh.quad_weights, "all");
            xMean = 1/Area*sum(regionMesh.geo_map.*permute(regionMesh.jacdet, [3,1,2]).*permute(regionMesh.quad_weights, [3,1,2]), [2,3]);
            if snapNearest
                [~, indMinDist] = min(vecnorm(obj.Nodes(regionNodeNumbers,:)-xMean', 2, 2));
                xCenter = regionNodes(indMinDist, 1:2);
            else
                xCenter = xMean';
            end
            xCenter = xCenter*[cos(obj.RotationAngle), sin(obj.RotationAngle); -sin(obj.RotationAngle), cos(obj.RotationAngle)];
        end

        function calculatePlotValues(obj)
            for iRegion = 1:obj.NumberRegions
                obj.Geometry(iRegion).PlotDimensions = size(obj.Geometry(iRegion).Triangles);
                obj.Geometry(iRegion).PlotMeshScalar = MESH_2Dtri(obj.Nodes, obj.Geometry(iRegion).Triangles, "qw", [1;1;1], "qn", {[0,0,1],[0,1,0]});
                obj.Geometry(iRegion).PlotSpaceScalar = SPACE_2Dtri(obj.Geometry(iRegion).PlotMeshScalar, obj.Geometry(iRegion).Triangles);
            end
        end
        
        function plotRegionNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14;
                args.Color = "black";
                args.RegionNumbers (1,:) double = 1:obj.NumberRegions;
                args.SnapNearest (1,1) logical = false;  % true center (can be outside region) or snap to nearest patch center
            end
            obj.checkRotationUpdate();
            xmidRegions = zeros(2, numel(args.RegionNumbers));
            for iRegion = 1:numel(args.RegionNumbers)
                xmidRegions(:, iRegion) = obj.getRegionCenter(args.RegionNumbers(iRegion), args.SnapNearest);
            end
            text(xmidRegions(1,:), xmidRegions(2,:), string(args.RegionNumbers), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize, "Color", args.Color);
        end

        function plotBoundaryNumbers(obj, args)
            arguments
                obj 
                args.FontSize = 14;
                args.Color = "black";
            end
            obj.checkRotationUpdate();
            pos = zeros(numel(obj.Boundaries), 2);
            for iBnd = 1:numel(obj.Boundaries)
                nodePos = obj.NodesPlot(obj.Boundaries(iBnd).dofs, 1:2);
                meanPos = mean(nodePos, 1);
                [~, indMinDist] = min(vecnorm(nodePos-meanPos, 2, 2));
                pos(iBnd, :) = nodePos(indMinDist, 1:2);
            end
             text(pos(:,1), pos(:,2), string(1:numel(obj.Boundaries)), "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", args.FontSize, "Color", args.Color);
        end

        function plotBoundaryConditions(obj, t, FigureName, ValueName, args)
            arguments
                obj 
                t = 0;
                FigureName = "Geometry";
                ValueName = FigureName;
                args.LineWidthBoundaries = 1.5;
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkRotationUpdate(t);
                        
            for iBC = 1:numel(obj.BoundaryConditions)
                checkBC = obj.BoundaryConditions(iBC).BoundaryCondition;
                if isprop(checkBC, "BoundaryNumbers")
                    col = checkBC.getPlotColor();
                    bndName = "BC"+iBC;
                    obj.PlotTool.initializeLines(FigureName, ValueName, bndName, vertcat(obj.Boundaries(checkBC.BoundaryNumbers).PlotDimensions), args.FocusFigure);
                    
                    for iBnd =  1:numel(checkBC.BoundaryNumbers)
                        bndNr = checkBC.BoundaryNumbers(iBnd);
                        xValues = [obj.NodesPlot(obj.Boundaries(bndNr).lines(:,1), 1), obj.NodesPlot(obj.Boundaries(bndNr).lines(:,2), 1)];
                        yValues = [obj.NodesPlot(obj.Boundaries(bndNr).lines(:,1), 2), obj.NodesPlot(obj.Boundaries(bndNr).lines(:,2), 2)];
                        obj.PlotTool.setLines(FigureName, ValueName, bndName, iBnd, "XData", xValues, "YData", yValues);
                    end
                    obj.PlotTool.updateLines(FigureName, ValueName, bndName, "LineWidth", args.LineWidthBoundaries, "EdgeColor", col);
                end
            end
        end

        function checkRotationUpdate(obj, t)
            if exist("t", "var") && isfield(obj.Solution, "RotationAngle") && ~isscalar(obj.Solution.RotationAngle)
                obj.RotationAngle = interp1(obj.Solution.Time, obj.Solution.RotationAngle, t);
            end
            if obj.RotationAngle ~= obj.PlotAngle
                obj.PlotAngle = obj.RotationAngle;
                obj.NodesPlot = obj.Nodes*[cos(obj.PlotAngle), sin(obj.PlotAngle); -sin(obj.PlotAngle), cos(obj.PlotAngle)];
            end
        end
    end
end
