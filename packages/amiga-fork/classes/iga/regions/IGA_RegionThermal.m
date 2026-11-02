classdef IGA_RegionThermal < IGA_Region

    properties

    end

    methods

        function obj = IGA_RegionThermal(subdivisions, degrees, length, args)
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


        function postprocess(obj, ts, ys, ysder)
            obj.Solution.Time = ts;
            obj.Solution.Temperature = ys;
            obj.Solution.TemperaturePrime = ysder;
            if ~isempty(obj.EvaluationShapeValues)
                obj.Solution.EvaluationT = reshape(sum(full(obj.EvaluationShapeValues{1}).*permute(ys, [3,1,2]), 2), [], numel(ts));
            end
        end

        function plotScalar(obj, value, t, args)
            arguments
                obj IGA_RegionThermal;
                value {mustBeMember(value, ["Temperature", "EnergyDensity", "Conductivity"])};
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
                args.TempOffset = 0;
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
                    tpatch = obj.Solution.Temperature(dofs);
                    % tpatchPrime = obj.Solution.TemperaturePrime(dofs);
                else
                    tpatch = interp1(obj.Solution.Time, obj.Solution.Temperature(dofs,:)', t)';
                    % tpatchPrime = interp1(obj.Solution.Time, obj.Solution.TemperaturePrime(dofs,:)', t)';
                end
                % evaluation always at unrotated geometry
                pnts = obj.getGeoPlotPoints(obj.Geometry(patch).nurbs, obj.Geometry(patch).PlotPointsParam, 0);
                xEval = reshape(pnts{1,1}, [], 1); yEval = reshape(pnts{2,1}, [], 1);
                % choose correct value
                switch value
                    case "Temperature"
                        PlotValue = postprocess_thermal(obj.Geometry(patch).PlotSpaceScalar, tpatch, "T") - args.TempOffset;
                    case "EnergyDensity"
                        Temperature = postprocess_thermal(obj.Geometry(patch).PlotSpaceScalar, tpatch, "T");
                        PlotValue = mat.getEnergyDensity(Temperature, xEval, yEval);
                    case "Conductivity"
                        Temperature = postprocess_thermal(obj.Geometry(patch).PlotSpaceScalar, tpatch, "T");
                        PlotValue = mat.getCondutivity(Temperature, xEval, yEval);
                end
                obj.PlotTool.setScalar(args.FigureName, args.ValueName, args.StoreName, patch, "XData", obj.Geometry(patch).PlotPoints{1}, ...
                    "YData", obj.Geometry(patch).PlotPoints{2}, "ZData", PlotValue, "CData", PlotValue, "flat", args.FlattenPlot, "absolute", args.AbsoluteValues, "scaling", args.Scaling);
            end

            obj.PlotTool.updateScalar(args.FigureName, args.ValueName, args.StoreName, "Shading", args.Shading, "EdgeColor", args.EdgeColor);
        end


        function plotTemperature(obj, t, args)
            arguments
                obj
                t = 0
                args.cLimits = [-Inf, Inf];
                args.Limits {mustBeMember(args.Limits, ["auto", "manual", "method", "mode", "padded", "tickaligned", "tight"])} = "tight";
            end

            obj.plotScalar("Temperature", t, "FlattenPlot", true, "TempOffset", 273.15);
            hold on
            obj.plotBoundary(t, "FigureName", "Temperature", "ValueName", "Temperature", "ColorInterfaces", [1,1,1]*0.3);

            c = colorbar();
            c.TickLabelInterpreter = 'latex';
            c.Label.Interpreter = 'latex';
            c.Label.String = 'Temperature ($^\circ C$)';
            c.FontSize = 18;
            c.Label.FontSize = 20;
            % set(c, "Visible", "off");
            clim(args.cLimits);
            set(gca,'DataAspectRatio', [1 1 1e3])
            xlim(args.Limits);
            ylim(args.Limits);
            drawnow();
        end

    end
end