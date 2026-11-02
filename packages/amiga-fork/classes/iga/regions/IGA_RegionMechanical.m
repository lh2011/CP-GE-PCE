classdef IGA_RegionMechanical < IGA_Region

    properties

    end

    methods

        function obj = IGA_RegionMechanical(subdivisions, degrees, length, args)
            arguments
                subdivisions (1, 2) double  = [5, 5];
                degrees (1, 2) double = [2, 2];
                length (1, 1) double = 1;
                args.DegreeQuadrature = degrees + [1, 1];
            end
            obj@IGA_Region(subdivisions, degrees, length, "DegreeQuadrature", args.DegreeQuadrature);
        end

        function generateSpaces(obj)
            progressbar([obj.Name ': Generating function spaces: '])

            % pre-evaluated spaces
            obj.SpacesEval  = cell (1, obj.NumberPatches);

            % temporary cell arrays
            spaces  = cell (1, obj.NumberPatches);
            for iPatch = 1:obj.NumberPatches
                [knots, ~] = kntrefine (obj.Geometry(iPatch).nurbs.knots, obj.SubdivisionsPatches-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
                scalar_spaces_mech = cell (obj.Meshes.msh_patch{iPatch}.rdim, 1);
                for idim = 1:obj.Meshes.msh_patch{iPatch}.rdim
                    scalar_spaces_mech{idim} = sp_bspline(knots, obj.FormFunctionDegree, obj.Meshes.msh_patch{iPatch});
                    % scalar_spaces_mech{idim} = sp_nurbs(obj.Geometry(iPatch).nurbs, obj.Meshes.msh_patch{iPatch});
                end
                spaces{iPatch} = sp_vector(scalar_spaces_mech, obj.Meshes.msh_patch{iPatch});

                obj.SpacesEval{iPatch}  = sp_precompute(spaces{iPatch}, obj.MeshesEval{iPatch}, 'value', true, 'gradient', true);
                progressbar(iPatch/obj.NumberPatches);
            end
            obj.Spaces  = sp_multipatch(spaces, obj.Meshes, obj.PatchInterfaces, obj.BoundaryInterfaces);

            % init values
            obj.PlotValuesNeedUpdate = true;

            progressbar('finished');
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Displacement = ys;
        end

        function plotScalar(obj, value, t, args)
            arguments
                obj IGA_RegionMechanical;
                value {mustBeMember(value, ["Displacement", "DisplacementX", "DisplacementY", "StressVonMises", "StressX", "StressY", "EnergyDensity"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Scalar";
                args.EdgeColor  = "none";
                args.Shading {mustBeMember(args.Shading, ["interp", "flat"])} = "interp";
                args.Patches = 1:obj.NumberPatches;
                args.FlattenPlot (1,1) logical = false;
                args.AbsoluteValues (1,1) logical = false
                args.Scaling = 1;
                args.FocusFigure (1,1) logical = true;
                args.DisplacementScaling = 0;
            end
            obj.checkPlotValues();
            obj.checkRotationUpdate(t);

            if isempty(obj.Solution.Time)
                error("Run a solver before plotting results")
            end

            obj.PlotTool.initializeScalar(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensions), args.FocusFigure);

            for patch = args.Patches
                mat = obj.Regions(obj.Geometry(patch).Region).Material;
                dofs = obj.Spaces.gnum{patch};
                if isscalar(obj.Solution.Time)
                    upatch = obj.Solution.Displacement(dofs);
                else
                    upatch = interp1(obj.Solution.Time, obj.Solution.Displacement(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                pnts = obj.getGeoPlotPoints(obj.Geometry(patch).nurbs, obj.Geometry(patch).PlotPointsParam, 0);
                xEval = reshape(pnts{1,1}, [], 1); yEval = reshape(pnts{2,1}, [], 1);

                if args.DisplacementScaling ~= 0
                    ux = args.DisplacementScaling*postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementX");
                    uy = args.DisplacementScaling*postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementY");
                else
                    ux = zeros(size(xEval));
                    uy = zeros(size(yEval));
                end
                [x, y] = rotdata(xEval+ux, yEval+uy, obj.PlotAngle);

                % choose correct value
                switch value
                    case "Displacement"
                        PlotValue = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "Displacement");
                    case "DisplacementX"
                        PlotValue = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementX");
                    case "DisplacementY"
                        PlotValue = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementY");
                    case "StressVonMises"
                        strain = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "Strain");
                        stress = mat.getStress(strain);
                        sigxx = stress(1,:);
                        sigyy = stress(2,:);
                        sigxy = stress(3,:);
                        PlotValue = sqrt(sigxx.*sigxx+sigyy.*sigyy-sigxx.*sigyy+3*sigxy.*sigxy) / 1e6;
                    case "StressX" 
                        strain = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "Strain");
                        stress = mat.getStress(strain);
                        PlotValue = stress(1,:) / 1e6;
                    case "StressY"
                        strain = postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "Strain");
                        stress = mat.getStress(strain);
                        PlotValue = stress(2,:) / 1e6;
                end
                obj.PlotTool.setScalar(args.FigureName, args.ValueName, args.StoreName, patch, "XData", x, ...
                    "YData", y, "ZData", PlotValue, "CData", PlotValue, "flat", args.FlattenPlot, "absolute", args.AbsoluteValues, "scaling", args.Scaling);
            end

            obj.PlotTool.updateScalar(args.FigureName, args.ValueName, args.StoreName, "Shading", args.Shading, "EdgeColor", args.EdgeColor);
        end

        function plotDisplacement(obj, t, args)
            arguments
                obj
                t = 0
                args.cLimits = [-Inf, Inf];
                args.Limits {mustBeMember(args.Limits, ["auto", "manual", "method", "mode", "padded", "tickaligned", "tight"])} = "tight";
            end
            % obj.plotScalar("Displacement", t, "FlattenPlot", true);
            % hold on
            % obj.plotBoundary(t, "FigureName", "Displacement", "ValueName", "Displacement", "ColorInterfaces", [1,1,1]*0.3);
            % 
            % c = colorbar();
            % c.TickLabelInterpreter = 'latex';
            % c.Label.Interpreter = 'latex';
            % c.Label.String = 'TBD ($^\circ C$)';
            % c.FontSize = 18;
            % c.Label.FontSize = 20;
            % % set(c, "Visible", "off");
            % clim(args.cLimits);
            % set(gca,'DataAspectRatio', [1 1 1e3])
            % xlim(args.Limits);
            % ylim(args.Limits);
            % drawnow();
        end

        % overwrite standard function for 2-dimensional spaces
        function calculatePlotValues(obj)
            progressbar([obj.Name ': Calculating function space at plot points: '])
            for iPatch = 1:obj.NumberPatches
                for iDim = 1:2
                    [obj.Geometry(iPatch).PlotMeshScalar, spaces_scalar{iDim}] = sp_eval_pnts(obj.Spaces.sp_patch{iPatch}.scalar_spaces{iDim}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsParam); 
                end
                obj.Geometry(iPatch).PlotSpaceScalar = sp_vector(spaces_scalar, obj.Geometry(iPatch).PlotMeshScalar);
                progressbar(iPatch/obj.NumberPatches);
            end
            % TBD boundaries and interfaces??
            progressbar('finished');
        end

        % overwrite standard function for 2-dimensional spaces
        function calculatePlotValuesVector(obj)
            progressbar([obj.Name ': Calculating function space at vector plot points: '])
            for iPatch = 1:obj.NumberPatches
                for iDim = 1:2
                    [obj.Geometry(iPatch).PlotMeshVector, spaces_scalar{iDim}] = sp_eval_pnts(obj.Spaces.sp_patch{iPatch}.scalar_spaces{iDim}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsVectorParam); 
                end
                obj.Geometry(iPatch).PlotSpaceVector = sp_vector(spaces_scalar, obj.Geometry(iPatch).PlotMeshVector);
                progressbar(iPatch/obj.NumberPatches);
            end
            progressbar('finished');
        end

        % Overwrite standard function for 2-dimensional spaces
        function [leftDOFsX, leftDOFsY, rightDOFsX, rightDOFsY] = getBoundaryDOFsperiodic(obj, boundaryNumbersLeft, boundaryNumbersRight)
            assert(numel(boundaryNumbersLeft) == numel(boundaryNumbersRight), "Number of left and right sides must be the same!");
            leftDOFsX = [];
            leftDOFsY = [];
            rightDOFsX = [];
            rightDOFsY = [];
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
                    leftDOFsX = union(leftDOFsX, globalDOFsLeft(1:end/2), 'stable');
                    leftDOFsY = union(leftDOFsY, globalDOFsLeft((end/2)+1:end), 'stable');
                    rightDOFsX = union(rightDOFsX, globalDOFsRight(1:end/2), 'stable');
                    rightDOFsY = union(rightDOFsY, globalDOFsRight((end/2)+1:end), 'stable');
                elseif abs(vecnorm(nrbeval(nrbRight, 0)) - vecnorm(nrbeval(nrbLeft, 1))) < tol
                    leftDOFsX = union(leftDOFsX, globalDOFsLeft(1:end/2), 'stable');
                    leftDOFsY = union(leftDOFsY, globalDOFsLeft((end/2)+1:end), 'stable');
                    rightDOFsX = union(rightDOFsX, flip(globalDOFsRight(1:end/2)), 'stable');
                    rightDOFsY = union(rightDOFsY, flip(globalDOFsRight((end/2)+1:end)), 'stable');
                else
                    error("Boundaries are not matching!")
                end
            end
        end
    
        % Overwrite standard function for deformation
        function plotGeometry(obj, t, args)
            arguments
                obj
                t = 0;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
                args.FocusFigure (1,1) logical = true;
                args.DisplacementScaling = 0;
            end
            
            obj.plotRegions(t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", args.StoreName, "FocusFigure", args.FocusFigure, "DisplacementScaling", args.DisplacementScaling);
            obj.plotBoundary(t, "FigureName", args.FigureName, "ValueName", args.ValueName, "FocusFigure", args.FocusFigure, "DisplacementScaling", args.DisplacementScaling);
            daspect([1,1,1]);
            % axis equal
        end


        % Overwrite standard function for deformation
        function plotRegions(obj, t, args)
            arguments
                obj
                t = 0;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
                args.FocusFigure (1,1) logical = true;
                args.DisplacementScaling = 0;
            end
            obj.checkRotationUpdate(t);

            obj.PlotTool.initializeColor(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensions), args.FocusFigure);

            for iRegion = 1:obj.NumberRegions
                patches = obj.Regions(iRegion).Patches;
                if isempty(obj.Regions(iRegion).Material)
                    col = "white";
                    alph = 1;
                    warning("IGA_RegionMechanical:MaterialNotDefined: Material for region %i is not yet defined!", iRegion)
                else
                    col = obj.Regions(iRegion).Material.getPlotColor();
                    alph = obj.Regions(iRegion).Material.getPlotAlpha();
                end
                for patch = patches
                    xPlot = reshape(obj.Geometry(patch).PlotPoints{1}, [], 1); yPlot = reshape(obj.Geometry(patch).PlotPoints{2}, [], 1);
                    if args.DisplacementScaling ~= 0
                        dofs = obj.Spaces.gnum{patch};
                        if isscalar(obj.Solution.Time)
                            upatch = obj.Solution.Displacement(dofs);
                        else
                            upatch = interp1(obj.Solution.Time, obj.Solution.Displacement(dofs,:)', t)';
                        end
                        ux = args.DisplacementScaling*postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementX");
                        uy = args.DisplacementScaling*postprocess_mechanical(obj.Geometry(patch).PlotSpaceScalar, upatch, "DisplacementY");
                        [ux, uy] = rotdata(ux, uy, obj.PlotAngle);
                    else
                        ux = zeros(size(xPlot));
                        uy = zeros(size(yPlot));
                    end
                    xValues = reshape(xPlot+ux, size(obj.Geometry(patch).PlotPoints{1}));
                    yValues = reshape(yPlot+uy, size(obj.Geometry(patch).PlotPoints{2}));

                    obj.PlotTool.setColor(args.FigureName, args.ValueName, args.StoreName, patch, "XData", xValues, ...
                        "YData", yValues, "Color", col, "Alpha", alph);
                end
            end 
            drawOrder = horzcat(obj.Regions.Patches);
            obj.PlotTool.updateColor(args.FigureName, args.ValueName, args.StoreName, drawOrder);
        end

        function [uxLine, uyLine] = getBoundaryDisplacement(obj, patchNr, side, scaling)
            if scaling == 0
                uxLine = 0;
                uyLine = 0;
                return
            end
            dofs = obj.Spaces.gnum{patchNr};
            if isscalar(obj.Solution.Time)
                upatch = obj.Solution.Displacement(dofs);
            else
                upatch = interp1(obj.Solution.Time, obj.Solution.Displacement(dofs,:)', t)';
            end
            ux = scaling*postprocess_mechanical(obj.Geometry(patchNr).PlotSpaceScalar, upatch, "DisplacementX");
            uy = scaling*postprocess_mechanical(obj.Geometry(patchNr).PlotSpaceScalar, upatch, "DisplacementY");
            [ux, uy] = rotdata(ux, uy, obj.PlotAngle);
            % patch values
            uxPatch = reshape(ux, obj.Geometry(patchNr).PlotDimensions);
            uyPatch = reshape(uy, obj.Geometry(patchNr).PlotDimensions);
            % extract side
            switch side
                case 1
                    xBnd = uxPatch(1, :); yBnd = uyPatch(1,:);
                case 2
                    xBnd = uxPatch(end, :); yBnd = uyPatch(end, :);
                case 3
                    xBnd = uxPatch(:, 1); yBnd = uyPatch(:, 1);
                case 4
                    xBnd = uxPatch(:, end); yBnd = uyPatch(:, end);
            end
            uxLine = [reshape(xBnd(1:end-1), [], 1), reshape(xBnd(2:end), [], 1)];
            uyLine = [reshape(yBnd(1:end-1), [], 1), reshape(yBnd(2:end), [], 1)];
        end

        % Overwrite standard function for deformation
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
                args.DisplacementScaling = 0;
            end
            obj.checkRotationUpdate(t);
            if ~isempty(obj.PatchInterfaces)
                % plot patch interfaces with no material jump
                InterfacesNoMat = obj.PatchInterfaces(vertcat(obj.PatchInterfaces.IsMaterialInterface) == 0);
                if ~isempty(InterfacesNoMat)
                    obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "InterfacesNoMat", vertcat(InterfacesNoMat.PlotDimensions), args.FocusFigure);
                    for iInt = 1:numel(InterfacesNoMat)
                        patchNr = InterfacesNoMat(iInt).patch1;
                        side = InterfacesNoMat(iInt).side1;
                        [ux, uy] = obj.getBoundaryDisplacement(patchNr, side, args.DisplacementScaling);
                        obj.PlotTool.setLines(args.FigureName, args.ValueName, "InterfacesNoMat", iInt, ...
                            "XData", InterfacesNoMat(iInt).PlotPointsX+ux, "YData", InterfacesNoMat(iInt).PlotPointsY+uy);
                    end
                    obj.PlotTool.updateLines(args.FigureName, args.ValueName, "InterfacesNoMat", "LineWidth", args.LineWidthInterfaces, "EdgeColor", args.ColorInterfaces);
                end
    
                % plot patch interfaces with material jump
                InterfacesMat = obj.PatchInterfaces(vertcat(obj.PatchInterfaces.IsMaterialInterface) == 1);
                if ~isempty(InterfacesMat)
                    obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "InterfacesMat", vertcat(InterfacesMat.PlotDimensions), args.FocusFigure);
                    for iInt = 1:numel(InterfacesMat)
                        patchNr = InterfacesMat(iInt).patch1;
                        side = InterfacesMat(iInt).side1;
                        [ux, uy] = obj.getBoundaryDisplacement(patchNr, side, args.DisplacementScaling);
                        obj.PlotTool.setLines(args.FigureName, args.ValueName, "InterfacesMat", iInt, ...
                            "XData", InterfacesMat(iInt).PlotPointsX+ux, "YData", InterfacesMat(iInt).PlotPointsY+uy);
                    end
                    obj.PlotTool.updateLines(args.FigureName, args.ValueName, "InterfacesMat", "LineWidth", args.LineWidthMaterials, "EdgeColor", args.ColorMaterials);
                end
            end
            % plot boundaries
            obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "Boundaries", vertcat(obj.PatchBoundaries.PlotDimensions), args.FocusFigure);
            for iBnd = 1:numel(obj.PatchBoundaries)
                patchNr = obj.PatchBoundaries(iBnd).patches;
                side = obj.PatchBoundaries(iBnd).faces;
                [ux, uy] = obj.getBoundaryDisplacement(patchNr, side, args.DisplacementScaling);
                obj.PlotTool.setLines(args.FigureName, args.ValueName, "Boundaries", iBnd, ...
                    "XData", obj.PatchBoundaries(iBnd).PlotPointsX+ux, "YData", obj.PatchBoundaries(iBnd).PlotPointsY+uy);
            end
            obj.PlotTool.updateLines(args.FigureName, args.ValueName, "Boundaries", "LineWidth", args.LineWidthBoundaries, "EdgeColor", args.ColorBoundaries);
        end
    
    end
end