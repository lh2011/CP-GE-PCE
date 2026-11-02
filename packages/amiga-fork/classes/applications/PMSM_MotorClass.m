classdef PMSM_MotorClass < EQS_System
    
    properties
        Rotor, Stator;
        Coupling;
    end
    
    methods
        function obj = PMSM_MotorClass(rotor, stator, coupling)
            obj.Rotor = rotor;
            obj.Stator = stator;
            obj.Coupling = coupling;
            obj.addElement(rotor);
            obj.addElement(stator);
            obj.addElement(coupling)
        end

        function setRotationAngle(obj, angle, unit)
            arguments
                obj 
                angle 
                unit {mustBeMember(unit, ["rad", "deg"])} = "rad";
            end
            if unit == "deg"
                angle = deg2rad(angle);
            end
            obj.Coupling.setRotationAngle(angle);
        end

        function plotGeometry(obj, t, plotRemanence, args)
            arguments
                obj
                t = 0;
                plotRemanence (1,1) logical = true;
                args.FigureName = "Geometry";
                args.ValueName = "Geometry";
                args.StoreName = "Geometry";
            end
            unpackedArgs = namedargs2cell(args);

            obj.Rotor.plotRegions(t, unpackedArgs{:});
            obj.Rotor.plotBoundary(t, "FigureName", args.FigureName, "ValueName", args.ValueName);
            obj.Stator.plotRegions(t, unpackedArgs{:});
            obj.Stator.plotBoundary(t, "FigureName", args.FigureName, "ValueName", args.ValueName);
            if plotRemanence
                obj.Rotor.plotVector("Remanence", t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", "Remanence", LineWidth=2);
                obj.Stator.plotVector("Remanence", t, "FigureName", args.FigureName, "ValueName", args.ValueName, "StoreName", "Remanence", LineWidth=2);
            end
            axis equal
        end

        function plotMagneticPotentialLines(obj, t, args)
            arguments
                obj 
                t = 0;
                args.nLevels = 20; 
            end
            if isscalar(obj.Solution.Time)
                allValues = [obj.Rotor.Solution.Values; obj.Stator.Solution.Values];
            else
                allValues = [interp1(obj.Rotor.Solution.Time, obj.Rotor.Solution.Values', t)'; ...
                    interp1(obj.Stator.Solution.Time, obj.Stator.Solution.Values', t)'];
            end
            limits = [min(allValues), max(allValues)];
            limits = [min([limits, -limits]), max([limits, -limits])];

            obj.Rotor.plotMagneticPotentialLines(t, "nLevels", args.nLevels, "LevelLimits", limits);
            obj.Stator.plotMagneticPotentialLines(t, "nLevels", args.nLevels, "LevelLimits", limits);
        end

        function plotMagneticFluxDensity(obj, t, args)
            arguments
                obj 
                t = 0
                args.cLimits = [0, 2];
                args.Scaling = 0.5;
                args.Color = "black";
                args.LineWidth = 0.5;
                args.Limits {mustBeMember(args.Limits, ["auto", "manual", "method", "mode", "padded", "tickaligned", "tight"])} = "tight";
            end
            
            obj.Rotor.plotScalar("MagneticFluxDensity", t, "FlattenPlot", true);
            hold on
            obj.Stator.plotScalar("MagneticFluxDensity", t, "FlattenPlot", true);
            
            obj.Rotor.plotBoundary(t, "FigureName", "MagneticFluxDensity", "ValueName", "MagneticFluxDensity");%, "ColorInterfaces", [1,1,1]*0.3);
            obj.Stator.plotBoundary(t, "FigureName", "MagneticFluxDensity", "ValueName", "MagneticFluxDensity");%, "ColorInterfaces", [1,1,1]*0.3);
            obj.Rotor.plotVector("MagneticFluxDensity", t, Scaling=args.Scaling, Color=args.Color, LineWidth=args.LineWidth);
            obj.Stator.plotVector("MagneticFluxDensity", t, Scaling=args.Scaling, Color=args.Color, LineWidth=args.LineWidth);
            
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
    end
end
