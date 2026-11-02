classdef FEM_RegionMagnetic < FEM_Region

    properties
        PlotValuesNeedUpdate = true;
    end

    methods
        function obj = FEM_RegionMagnetic(length, args)
            arguments
                length (1, 1) double = 1;
                args.degreeQuadrature = 1;
            end
            obj@FEM_Region(length, "degreeQuadrature", args.degreeQuadrature);
        end

        function generateSpaces(obj)
            mytimer("Spaces", "start");
            for iRegion = 1:numel(obj.Regions)
                obj.SpacesEval{iRegion} = SPACE_2Dtri(obj.MeshesEval{iRegion}, obj.Geometry(iRegion).Triangles);
                obj.Spaces.gnum{iRegion} = obj.SpacesEval{iRegion}.gnum;
                obj.Spaces.ncomp = 1;
            end
            obj.Spaces.ndof = size(obj.Nodes, 1);
            obj.NumberDOF = obj.Spaces.ndof;
            
            mytimer("Spaces", "stop");
        end

        function plotGeometry(obj, t, plotRemanence, args)
            arguments
                obj
                t = 0;
                plotRemanence (1,1) logical = true;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
                args.EdgeColor = "none";
                args.LineWidth = 0.5;
            end
            args1 = namedargs2cell(args);
            plotGeometry@FEM_Region(obj, t, args1{:});
            if plotRemanence
                obj.plotVector("Remanence", t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", "Remanence", LineWidth=2);
            end
        end

        function plotMagneticFluxDensity(obj, t, args)
            arguments
                obj
                t = 0
                args.cLimits = [0, 2];
                args.Scaling = 0.5;
                args.ArrowColor = "black";
                args.LineWidth = 0.5;
                args.EdgeColor = "none";
                args.Limits {mustBeMember(args.Limits, ["auto", "manual", "method", "mode", "padded", "tickaligned", "tight"])} = "tight";
            end

            obj.plotScalar("MagneticFluxDensity", t, "FlattenPlot", true, "EdgeColor", args.EdgeColor);
            % obj.plotContour("MagneticFluxDensity", t, "filled", true, "LineColor", "none", "LevelLimits", args.cLimits);
            hold on
            obj.plotBoundary(t, "MagneticFluxDensity", "MagneticFluxDensity");
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
        end

        function plotMagneticPotentialLines(obj, t, args)
            arguments
                obj;
                t = 0;
                args.nLevels = 20;
                args.LevelLimits = [];
            end
            obj.plotGeometry(t, false, "FigureName", "MagneticPotentialLines", "ValueName", "MagneticPotentialLines", "StoreName", "Geometry");
            obj.plotContour("MagneticVectorPotential", t, "FigureName", "MagneticPotentialLines", "nLevels", args.nLevels, "LevelLimits", args.LevelLimits);
        end
        
        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Values = ys;
            obj.Solution.ValuesDer = ysprime;
            % point evaluations
            % if ~isempty(obj.EvaluationShapeValues)
            %     obj.EvaluationValues = squeeze(sum((sum(double(obj.EvaluationShapeValues{2}).*permute(ys, [4,3,1,2]), 3)).^2, 2).^0.5);
            %     obj.Solution.EvaluationAz = obj.EvaluationShapeValues{1}*ys;
            % end
        end


        function plotScalar(obj, value, t, args)
            arguments
                obj FEM_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticVectorPotential", ...
                    "MagneticFieldStrength", "Remanence", "Permeability", "Reluctivity", "EnergyDensity", ...
                    "LossDensity", "CurrentDensity", "Demagnetization"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Scalar";
                args.EdgeColor  = "none";
                args.Shading {mustBeMember(args.Shading, ["interp", "flat"])} = "interp";
                args.Regions = 1:obj.NumberRegions;
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
            
            for region = args.Regions
                mat = obj.Regions(region).Material;
                dofs = obj.Spaces.gnum{region};
                if isscalar(obj.Solution.Time)
                    apatch = obj.Solution.Values(dofs);
                    apatchPrime = obj.Solution.ValuesDer(dofs);
                else
                    apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    apatchPrime = interp1(obj.Solution.Time, obj.Solution.ValuesDer(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                xEval = reshape(obj.Geometry(region).PlotMeshScalar.geo_map(1,:,:), [], 1);
                yEval = reshape(obj.Geometry(region).PlotMeshScalar.geo_map(2,:,:), [], 1);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValue = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                    case "MagneticVectorPotential"
                        PlotValue = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Az");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        [Hx, Hy] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                        PlotValue = (Hx.^2 + Hy.^2).^0.5;
                    case "Remanence"
                        [Brx, Bry] = mat.getRemanence(xEval, yEval);
                        PlotValue = (Brx.^2 + Bry.^2).^0.5;
                    case "Permeability"
                        Bmag = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getPermeability(Bmag)/(4*pi*1e-7);
                    case "Reluctivity"
                        Bmag = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getReluctivity(Bmag)*(4*pi*1e-7);
                    case "EnergyDensity"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getEnergyDensity(Bx, By, xEval, yEval);
                    case "CurrentDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getCurrentDensity(t, AzPrime, region, xEval, yEval);
                    case "LossDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getLossDensity(t, AzPrime, region, xEval, yEval);
                    case "Demagnetization"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getDemagnetization(Bx, By, xEval, yEval);
                end
                xRot = [xEval, yEval]*[cos(obj.PlotAngle), sin(obj.PlotAngle); -sin(obj.PlotAngle), cos(obj.PlotAngle)];
                obj.PlotTool.setScalar(args.FigureName, args.ValueName, args.StoreName, region, "XData", xRot(:,1), ...
                    "YData", xRot(:,2), "ZData", PlotValue, "CData", PlotValue, "flat", args.FlattenPlot, "absolute", args.AbsoluteValues, "scaling", args.Scaling);
            end

            obj.PlotTool.updateScalar(args.FigureName, args.ValueName, args.StoreName, "Shading", args.Shading, "EdgeColor", args.EdgeColor);
        end

        function plotVector(obj, value, t, args)
            arguments
                obj FEM_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticFieldStrength", "Remanence"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Vector";
                args.Scaling = 0.5;
                args.Color = "black";
                args.LineWidth = 0.5;
                args.Regions = 1:obj.NumberRegions;
                args.FocusFigure (1,1) logical = true;
            end

            obj.checkPlotValues();
            obj.checkRotationUpdate(t);
            
            if ~strcmp(value, "Remanence")
                if isempty(obj.Solution.Time)
                    error("Run a solver before plotting results")
                end
            end

            obj.PlotTool.initializeVector(args.FigureName, args.ValueName, args.StoreName, vertcat(obj.Geometry.PlotDimensionsVector), args.FocusFigure);

            for region = args.Regions
                mat = obj.Regions(region).Material;
                if isempty(mat)
                    continue;
                end
                if ~strcmp(value, "Remanence")
                    dofs = obj.Spaces.gnum{region};
                    if isscalar(obj.Solution.Time)
                        apatch = obj.Solution.Values(dofs);
                    else
                        apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    end
                end
                % evaluation always at unrotated geometry
                xEval = reshape(obj.Geometry(region).PlotMeshVector.geo_map(1,:,:), [], 1);
                yEval = reshape(obj.Geometry(region).PlotMeshVector.geo_map(2,:,:), [], 1);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValueU = postprocess_magnetic(obj.Geometry(region).PlotSpaceVector, apatch, "Bx");
                        PlotValueV = postprocess_magnetic(obj.Geometry(region).PlotSpaceVector, apatch, "By");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceVector, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceVector, apatch, "By");
                        [PlotValueU, PlotValueV] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                    case "Remanence"
                        [PlotValueU, PlotValueV] = mat.getRemanence(xEval, yEval);
                end
                xRot = [xEval, yEval]*[cos(obj.PlotAngle), sin(obj.PlotAngle); -sin(obj.PlotAngle), cos(obj.PlotAngle)];
                obj.PlotTool.setVector(args.FigureName, args.ValueName, args.StoreName, region, "XData", xRot(:,1), ...
                    "YData", xRot(:,2), "UData", PlotValueU, "VData", PlotValueV, "Rotation", obj.RotationAngle);
            end

            obj.PlotTool.updateVector(args.FigureName, args.ValueName, args.StoreName, ...
                Scaling=args.Scaling, Color=args.Color, LineWidth=args.LineWidth);
        end

        function plotContour(obj, value, t, args)
            arguments
                obj FEM_RegionMagnetic;
                value {mustBeMember(value, ["MagneticFluxDensity", "MagneticVectorPotential", ...
                    "MagneticFieldStrength", "Remanence", "Permeability", "Reluctivity", "EnergyDensity", ...
                    "LossDensity", "CurrentDensity", "Demagnetization"])};
                t=0;
                args.FigureName = value;
                args.ValueName = value;
                args.StoreName = "Contour";
                args.LevelLimits = [];
                args.nLevels = 20;
                args.LineColor = "black";
                args.LineWidth = 0.5;
                % args.Filled (1,1) logical = true;
                args.Regions = 1:obj.NumberRegions;
                args.AbsoluteValues (1,1) logical = false;
                args.Scaling = 1;
                args.FocusFigure (1,1) logical = true;
            end
            obj.checkPlotValues();
            obj.checkRotationUpdate(t);

            if isempty(obj.Solution.Time)
                error("Run a solver before plotting results")
            end
            if isempty(args.LevelLimits)
                if isscalar(obj.Solution.Time)
                    args.LevelLimits = [min(obj.Solution.Values), max(obj.Solution.Values)];
                else
                    a = interp1(obj.Solution.Time, obj.Solution.Values(:,:)', t)';
                    args.LevelLimits = [min(a), max(a)];
                end
            end
            levels = linspace(args.LevelLimits(1), args.LevelLimits(2), args.nLevels);
            % Contour is set to be made of three lines per triangle
            nLines = vertcat(obj.Geometry.PlotDimensions)*args.nLevels; nLines(:,2) = 2;
           
            obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "ContourLines1", nLines, args.FocusFigure);
            obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "ContourLines2", nLines, args.FocusFigure);
            obj.PlotTool.initializeLines(args.FigureName, args.ValueName, "ContourLines3", nLines, args.FocusFigure);

            for region = args.Regions
                mat = obj.Regions(region).Material;
                dofs = obj.Spaces.gnum{region};
                if isscalar(obj.Solution.Time)
                    apatch = obj.Solution.Values(dofs);
                    apatchPrime = obj.Solution.ValuesDer(dofs);
                else
                    apatch = interp1(obj.Solution.Time, obj.Solution.Values(dofs,:)', t)';
                    apatchPrime = interp1(obj.Solution.Time, obj.Solution.ValuesDer(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                xEval = permute(obj.Geometry(region).PlotMeshScalar.geo_map(1,:,:), [2,3,1]);
                yEval = permute(obj.Geometry(region).PlotMeshScalar.geo_map(2,:,:), [2,3,1]);
                % choose correct value
                switch value
                    case "MagneticFluxDensity"
                        PlotValue = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                    case "MagneticVectorPotential"
                        PlotValue = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Az");
                    case "MagneticFieldStrength"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        [Hx, Hy] = mat.getMagneticFieldStrength(Bx, By, xEval, yEval);
                        PlotValue = (Hx.^2 + Hy.^2).^0.5;
                    case "Remanence"
                        [Brx, Bry] = mat.getRemanence(xEval, yEval);
                        PlotValue = (Brx.^2 + Bry.^2).^0.5;
                    case "Permeability"
                        Bmag = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getPermeability(Bmag)/(4*pi*1e-7);
                    case "Reluctivity"
                        Bmag = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bmag");
                        PlotValue = mat.getReluctivity(Bmag)*(4*pi*1e-7);
                    case "EnergyDensity"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getEnergyDensity(Bx, By, xEval, yEval);
                    case "CurrentDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getCurrentDensity(t, AzPrime, region, xEval, yEval);
                    case "LossDensity"
                        AzPrime = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatchPrime, "Az");
                        PlotValue = mat.getLossDensity(t, AzPrime, region, xEval, yEval);
                    case "Demagnetization"
                        Bx = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "Bx");
                        By = postprocess_magnetic(obj.Geometry(region).PlotSpaceScalar, apatch, "By");
                        PlotValue = mat.getDemagnetization(Bx, By, xEval, yEval);
                end
                
                [xEval, yEval] = rotdata(xEval, yEval, obj.PlotAngle);

                zEval = reshape(PlotValue, 3, [])';
                v1 = [xEval(1,:)', yEval(1,:)', zEval(:,1)];
                v2 = [xEval(2,:)', yEval(2,:)', zEval(:,2)];
                v3 = [xEval(3,:)', yEval(3,:)', zEval(:,3)];
                secPoints = calcTriangleIntersections(v1, v2, v3, levels);

                xdata = [reshape(secPoints(1,:,:,1), [], 1), reshape(secPoints(1,:,:,2), [], 1)];
                ydata = [reshape(secPoints(2,:,:,1), [], 1), reshape(secPoints(2,:,:,2), [], 1)];
                obj.PlotTool.setLines(args.FigureName, args.ValueName,"ContourLines1", region,"XData", xdata, "YData", ydata);
                xdata = [reshape(secPoints(1,:,:,2), [], 1), reshape(secPoints(1,:,:,3), [], 1)];
                ydata = [reshape(secPoints(2,:,:,2), [], 1), reshape(secPoints(2,:,:,3), [], 1)];
                obj.PlotTool.setLines(args.FigureName, args.ValueName,"ContourLines2", region,"XData", xdata, "YData", ydata);
                xdata = [reshape(secPoints(1,:,:,3), [], 1), reshape(secPoints(1,:,:,1), [], 1)];
                ydata = [reshape(secPoints(2,:,:,3), [], 1), reshape(secPoints(2,:,:,1), [], 1)];
                obj.PlotTool.setLines(args.FigureName, args.ValueName,"ContourLines3", region,"XData", xdata, "YData", ydata);
            end
            obj.PlotTool.updateLines(args.FigureName, args.ValueName,"ContourLines1", "LineWidth", args.LineWidth, "EdgeColor", args.LineColor);
            obj.PlotTool.updateLines(args.FigureName, args.ValueName,"ContourLines2", "LineWidth", args.LineWidth, "EdgeColor", args.LineColor);
            obj.PlotTool.updateLines(args.FigureName, args.ValueName,"ContourLines3", "LineWidth", args.LineWidth, "EdgeColor", args.LineColor);
        end

        
        function calculatePlotValuesVector(obj)
            for iRegion = 1:obj.NumberRegions
                obj.Geometry(iRegion).PlotDimensionsVector = [size(obj.Geometry(iRegion).Triangles, 1), 1];
                obj.Geometry(iRegion).PlotMeshVector = MESH_2Dtri(obj.Nodes, obj.Geometry(iRegion).Triangles, "qw", 1, "qn", {0.5, 0.5});
                obj.Geometry(iRegion).PlotSpaceVector = SPACE_2Dtri(obj.Geometry(iRegion).PlotMeshVector, obj.Geometry(iRegion).Triangles);
            end
        end

        function checkPlotValues(obj)
            if obj.PlotValuesNeedUpdate
                obj.calculatePlotValues();
                obj.calculatePlotValuesVector();
            end
            obj.PlotValuesNeedUpdate = false;
        end

    end
end