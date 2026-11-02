% class for holding the X vector
% for current density excitation and field circuit coupling

classdef MNA_CircuitInterface < handle
    properties (Access = public)
        ParentRegion;           % Region reference
        ParentMaterial;         % Material reference
        Phases;                 % Definition for phases, corresponding patches and direction
        NumberWindings;         % Number of turns
        A = [1; -1];            % Incidence vector for circuit coupling
        Sigma;                  % Conductivity of conductor
        PlotColor;              % Color of winding plot
    end

    methods (Access = public)
        function obj = MNA_CircuitInterface(material, phaseDefinition, args)
            arguments
                material 
                phaseDefinition 
                args.NumberWindings = 1;
                args.PlotColor = [0, 0, 0];
            end
            obj.PlotColor = args.PlotColor;
            obj.ParentMaterial = material;
            obj.ParentRegion = material.Region;
            obj.NumberWindings = args.NumberWindings;
            obj.Sigma = material.Sigma;
            for iPhase = 1:numel(phaseDefinition)
                obj.Phases(iPhase).Regions = phaseDefinition(iPhase).Regions;
                obj.Phases(iPhase).Direction = phaseDefinition(iPhase).Direction;
                mesh = material.extractMeshSpace(obj.Phases(iPhase).Regions);
                obj.Phases(iPhase).Area = sum(mesh.jacdet.*mesh.quad_weights, "all");
                obj.Phases(iPhase).Conductance = obj.Sigma*obj.Phases(iPhase).Area/obj.ParentRegion.Length;
            end
        end

        function X = getX(obj)
            X = zeros(obj.ParentMaterial.Space.ndof, 1);
            for iPhase = 1:numel(obj.Phases)
                phaseRegions = obj.Phases(iPhase).Regions;
                [submesh, subspace, subindices] = obj.ParentMaterial.extractMeshSpace(phaseRegions);
                p = obj.Phases(iPhase).Direction;
                X(subindices) = X(subindices) + p*obj.NumberWindings/obj.Phases(iPhase).Area*submesh.length*int_f_v(subspace, submesh);
            end
        end

        function Q = getQ(obj)
            Q = zeros(obj.ParentMaterial.Space.ndof, 1);
            if numel(obj.Phases) > 1
                error("TBD, a single conductor that goes through the region has yet to be implemented");
            end
            for iPhase = 1:numel(obj.Phases)
                phaseRegions = obj.Phases(iPhase).Regions;
                [submesh, subspace, subindices] = obj.ParentMaterial.extractMeshSpace(phaseRegions);
                p = obj.Phases(iPhase).Direction;
                Q(subindices) = Q(subindices) + p*obj.Sigma*int_f_v(subspace, submesh);
            end
        end

        function Msigma = getMsigma(obj)
            Msigma = sparse(obj.ParentMaterial.Space.ndof, obj.ParentMaterial.Space.ndof);
            for iPhase = 1:numel(obj.Phases)
                phaseRegions = obj.Phases(iPhase).Regions;
                [submesh, subspace, subindices] = obj.ParentMaterial.extractMeshSpace(phaseRegions);
                Msigma(subindices, subindices) = Msigma(subindices, subindices) + submesh.length*obj.Sigma*int_u_v(subspace, subspace, submesh);
            end
        end
        
        function plotWindingCurrent(obj, args)
            arguments
                obj 
                args.RadialPlot (1,1) logical = true;   % rotate with position
                args.PlotColor = obj.PlotColor;
                args.PlotAlpha = 1;
                args.LineWidth = 3;
                args.ShrinkScale = 3;
                args.SnapNearest (1,1) logical = false;
            end
            
            for iPhase = 1:numel(obj.Phases)
                Xmean = obj.ParentRegion.getRegionCenter(horzcat(obj.Phases(iPhase).Regions), args.SnapNearest);
                if args.RadialPlot
                    angle = atan2(Xmean(2), Xmean(1));
                else
                    angle = 0;
                end
                if obj.Phases(iPhase).Direction == 1
                    plotInCurrent(Xmean(1), Xmean(2), (obj.Phases(iPhase).Area/pi)^0.5/args.ShrinkScale, args.PlotColor, args.PlotAlpha, args.LineWidth);
                else
                    plotOutCurrent(Xmean(1), Xmean(2), (obj.Phases(iPhase).Area/pi)^0.5/args.ShrinkScale, angle, args.PlotColor, args.PlotAlpha, args.LineWidth);
                end
            end

            function plotInCurrent(x, y, r, plotColor, alpha, lineWidth)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                p1 = nsidedpoly(1000, 'Center', [x y], 'Radius', r/5);
                plot(p, 'FaceColor', 'white', 'EdgeColor', plotColor, 'FaceAlpha', alpha, 'LineWidth', lineWidth)
                plot(p1, 'FaceColor', plotColor, 'EdgeColor', plotColor, 'FaceAlpha', alpha)
            end

            function plotOutCurrent(x, y, r, rot, plotColor, alpha, lineWidth)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                pts = [-r, -r, r, r; -r, r -r, r]/(2^0.5);
                R = [cos(rot), -sin(rot); sin(rot), cos(rot)];
                pts = R*pts;
                plot(p, 'FaceColor', 'white', 'EdgeColor', plotColor, 'FaceAlpha', alpha, 'LineWidth', lineWidth)
                plot(x + [pts(1, 1), pts(1, 4) ], y + [pts(2, 1), pts(2, 4) ], 'LineWidth', lineWidth, 'Color', plotColor)
                plot(x + [pts(1, 2), pts(1, 3) ], y + [pts(2, 2), pts(2, 3) ], 'LineWidth', lineWidth, 'Color', plotColor)
            end
        end
    end
end