classdef IGA_RegionMagnetic < IGA_Region

    properties
        
    end

    methods
        function obj = IGA_RegionMagnetic(subdivisions, degrees, length, args)
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
                spaces{iPatch}  = sp_bspline (knots, obj.FormFunctionDegree, obj.Meshes.msh_patch{iPatch});
                obj.SpacesEval{iPatch}  = sp_precompute(spaces{iPatch}, obj.MeshesEval{iPatch}, 'value', true, 'gradient', true);
                progressbar(iPatch/obj.NumberPatches);
            end
            obj.Spaces  = sp_multipatch(spaces, obj.Meshes, obj.PatchInterfaces, obj.BoundaryInterfaces);

            % init values
            obj.PlotValuesNeedUpdate = true;

            progressbar('finished');
        end

        function plotScalar(obj, value, t, args)
            arguments
                obj IGA_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticVectorPotential", ...
                    "MagneticFieldStrength", "Remanence", "Permeability", "Reluctivity", "EnergyDensity", ...
                    "LossDensity", "CurrentDensity", "Demagnetization"])};
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
            end
            obj.checkPlotValues();
            obj.checkRotationUpdate(t);

            if isempty(obj.Solution.Time)
                error("Run a solver before plotting results")
            end

            obj.PlotTool.initializeScalar(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensions), args.FocusFigure);

            for patch = args.Patches
                mat = obj.Regions(obj.Geometry(patch).Region).Material;
                reg = obj.Geometry(patch).Region;
                dofs = obj.Spaces.gnum{patch};
                if isscalar(obj.Solution.Time)
                    apatch = obj.Solution.Values(dofs);
                    apatchPrime = obj.Solution.ValuesDer(dofs);
                else
                    apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    apatchPrime = interp1(obj.Solution.Time, obj.Solution.ValuesDer(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                pnts = obj.getGeoPlotPoints(obj.Geometry(patch).nurbs, obj.Geometry(patch).PlotPointsParam, 0);
                xEval = reshape(pnts{1,1}, [], 1); yEval = reshape(pnts{2,1}, [], 1);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValue = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                    case "MagneticVectorPotential"
                        PlotValue = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Az");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        [Hx, Hy] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                        PlotValue = (Hx.^2 + Hy.^2).^0.5;
                    case "Remanence"
                        [Brx, Bry] = mat.getRemanence(xEval, yEval);
                        PlotValue = (Brx.^2 + Bry.^2).^0.5;
                    case "Permeability"
                        Bmag = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getPermeability(Bmag)/(4*pi*1e-7);
                    case "Reluctivity"
                        Bmag = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getReluctivity(Bmag)*(4*pi*1e-7);
                    case "EnergyDensity"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getEnergyDensity(Bx, By, xEval, yEval);
                    case "CurrentDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getCurrentDensity(t, AzPrime, reg, xEval, yEval);
                    case "LossDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getLossDensity(t, AzPrime, reg, xEval, yEval);
                    case "Demagnetization"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getDemagnetization(Bx, By, xEval, yEval);
                end
                obj.PlotTool.setScalar(args.FigureName, args.ValueName, args.StoreName, patch, "XData", obj.Geometry(patch).PlotPoints{1}, ...
                    "YData", obj.Geometry(patch).PlotPoints{2}, "ZData", PlotValue, "CData", PlotValue, "flat", args.FlattenPlot, "absolute", args.AbsoluteValues, "scaling", args.Scaling);
            end

            obj.PlotTool.updateScalar(args.FigureName, args.ValueName, args.StoreName, "Shading", args.Shading, "EdgeColor", args.EdgeColor);
        end

        function plotContour(obj, value, t, args)
            arguments
                obj IGA_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticVectorPotential", ...
                    "MagneticFieldStrength", "Remanence", "Permeability", "Reluctivity", "EnergyDensity", ...
                    "LossDensity", "CurrentDensity", "Demagnetization", "MagneticForceDensity"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Contour";
                args.LevelLimits = [];
                args.nLevels = 20;
                args.LineColor = "black";
                args.LineWidth = 0.5;
                args.Filled (1,1) logical = true;
                args.Patches = 1:obj.NumberPatches;
                args.AbsoluteValues (1,1) logical = false;
                args.Scaling = 1;
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkPlotValues();
            obj.checkRotationUpdate(t);

            if isempty(obj.Solution.Time)
                error("Run a solver before plotting results")
            end

            obj.PlotTool.initializeContour(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensions), args.FocusFigure);

            for patch = args.Patches
                mat = obj.Regions(obj.Geometry(patch).Region).Material;
                dofs = obj.Spaces.gnum{patch};
                if isscalar(obj.Solution.Time)
                    apatch = obj.Solution.Values(dofs);
                    apatchPrime = obj.Solution.ValuesDer(dofs);
                else
                    apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    apatchPrime = interp1(obj.Solution.Time, obj.Solution.ValuesDer(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                pnts = obj.getGeoPlotPoints(obj.Geometry(patch).nurbs, obj.Geometry(patch).PlotPointsParam, 0);
                xEval = reshape(pnts{1,1}, [], 1); yEval = reshape(pnts{2,1}, [], 1);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValue = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                    case "MagneticVectorPotential"
                        PlotValue = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Az");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        [Hx, Hy] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                        PlotValue = (Hx.^2 + Hy.^2).^0.5;
                    case "Remanence"
                        [Brx, Bry] = mat.getRemanence(xEval, yEval);
                        PlotValue = (Brx.^2 + Bry.^2).^0.5;
                    case "Permeability"
                        Bmag = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getPermeability(Bmag)/(4*pi*1e-7);
                    case "Reluctivity"
                        Bmag = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getReluctivity(Bmag)*(4*pi*1e-7);
                    case "EnergyDensity"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getEnergyDensity(Bx, By, xEval, yEval);
                    case "CurrentDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getCurrentDensity(t, AzPrime, patch, xEval, yEval);
                    case "LossDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getLossDensity(t, AzPrime, patch, xEval, yEval);
                    case "Demagnetization"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getDemagnetization(Bx, By, xEval, yEval);
                end
                obj.PlotTool.setContour(args.FigureName, args.ValueName, args.StoreName, patch, "XData", obj.Geometry(patch).PlotPoints{1}, ...
                    "YData", obj.Geometry(patch).PlotPoints{2}, "ZData", reshape(PlotValue, obj.Geometry(patch).PlotDimensions), "absolute", args.AbsoluteValues, "scaling", args.Scaling);
            end

            obj.PlotTool.updateContour(args.FigureName, args.ValueName, args.StoreName, ...
                LevelLimits=args.LevelLimits, nLevels=args.nLevels, LineColor=args.LineColor, LineWidth=args.LineWidth, filled=args.Filled);
        end

        function plotGeometry(obj, t, plotRemanence, args)
            arguments
                obj
                t = 0;
                plotRemanence (1,1) logical = true;
                args.FocusFigure (1,1) logical = true;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
            end
            args1 = namedargs2cell(args);
            plotGeometry@IGA_Region(obj, t, args1{:});
            if plotRemanence
                obj.plotVector("Remanence", t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", "Remanence", "LineWidth", 2, "FocusFigure", args.FocusFigure);
            end
        end

        function plotVector(obj, value, t, args)
            arguments
                obj IGA_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticFieldStrength", "Remanence"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Vector";
                args.Scaling = 0.5;
                args.Color = "black";
                args.LineWidth = 0.5;
                args.Patches = 1:obj.NumberPatches;
                args.FocusFigure (1,1) logical = true;
            end

            obj.checkRotationUpdate(t);
            if ~strcmp(value, "Remanence")
                obj.checkPlotValues();
                if isempty(obj.Solution.Time)
                    error("Run a solver before plotting results")
                end
            end

            obj.PlotTool.initializeVector(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensionsVector), args.FocusFigure);
            
            for patch = args.Patches
                mat = obj.Regions(obj.Geometry(patch).Region).Material;
                if isempty(mat)
                    continue;
                end
                if ~strcmp(value, "Remanence")
                    dofs = obj.Spaces.gnum{patch};
                    if isscalar(obj.Solution.Time)
                        apatch = obj.Solution.Values(dofs);
                    else
                        apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    end
                end
                % evaluation always at unrotated geometry
                pnts = obj.getGeoPlotPoints(obj.Geometry(patch).nurbs, obj.Geometry(patch).PlotPointsVectorParam, 0);
                xEval = reshape(pnts{1,1}, [], 1); yEval = reshape(pnts{2,1}, [], 1);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValueU = postprocess_magnetic(obj.Geometry(patch).PlotSpaceVector, apatch, "Bx");
                        PlotValueV = postprocess_magnetic(obj.Geometry(patch).PlotSpaceVector, apatch, "By");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(patch).PlotSpaceVector, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(patch).PlotSpaceVector, apatch, "By");
                        [PlotValueU, PlotValueV] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                    case "Remanence"
                        [PlotValueU, PlotValueV] = mat.getRemanence(xEval, yEval);
                end
                obj.PlotTool.setVector(args.FigureName, args.ValueName, args.StoreName, patch, "XData", obj.Geometry(patch).PlotPointsVector{1}, ...
                    "YData", obj.Geometry(patch).PlotPointsVector{2}, "UData", PlotValueU, "VData", PlotValueV, "Rotation", obj.RotationAngle);
            end

            obj.PlotTool.updateVector(args.FigureName, args.ValueName, args.StoreName, ...
                Scaling=args.Scaling, Color=args.Color, LineWidth=args.LineWidth);
        end

        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Values = ys;
            obj.Solution.ValuesDer = ysprime;
            % point evaluations
            if ~isempty(obj.EvaluationShapeValues)
                obj.Solution.EvaluationBmag = squeeze(sum((sum(double(obj.EvaluationShapeValues{2}).*permute(ys, [4,3,1,2]), 3)).^2, 2).^0.5);
                obj.Solution.EvaluationAz = obj.EvaluationShapeValues{1}*ys;
            end
        end

        %%%%%%%%%%%% PLOTTING

        function plotMagneticFluxDensity(obj, t, args)
            arguments
                obj
                t = 0
                args.cLimits = [0, 2];
                args.Scaling = 0.5;
                args.ArrowColor = "black";
                args.LineWidth = 0.5;
                args.Limits {mustBeMember(args.Limits, ["auto", "manual", "method", "mode", "padded", "tickaligned", "tight"])} = "tight";
            end

            obj.plotScalar("MagneticFluxDensity", t, "FlattenPlot", true);
            % obj.plotContour("MagneticFluxDensity", t, "filled", true, "LineColor", "none", "LevelLimits", args.cLimits);
            hold on
            obj.plotBoundary(t, "FigureName", "MagneticFluxDensity", "ValueName", "MagneticFluxDensity", "ColorInterfaces", [1,1,1]*0.3);
            obj.plotVector("MagneticFluxDensity", t, Scaling=args.Scaling, Color=args.ArrowColor, LineWidth=args.LineWidth);

            c = colorbar();
            c.TickLabelInterpreter = 'latex';
            c.Label.Interpreter = 'latex';
            c.Label.String = 'Magnetic flux density (T)';
            c.FontSize = 18;
            c.Label.FontSize = 20;
            % set(c, "Visible", "off");
            clim(args.cLimits);
            set(gca,'DataAspectRatio', [1 1 1e3])
            xlim(args.Limits);
            ylim(args.Limits);
            drawnow();
        end

        function plotMagneticPotentialLines(obj, t, args)
            arguments
                obj;
                t = 0;
                args.nLevels = 20;
                args.LevelLimits = [];
            end
            obj.plotGeometry(t, "FigureName", "MagneticPotentialLines", "ValueName", "MagneticPotentialLines", "StoreName", "Geometry");
            obj.plotContour("MagneticVectorPotential", t, "FigureName", "MagneticPotentialLines", "filled", false, "nLevels", args.nLevels, "LevelLimits", args.LevelLimits);
        end
    end
end