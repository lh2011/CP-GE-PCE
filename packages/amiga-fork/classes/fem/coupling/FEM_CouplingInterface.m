% class holding the meshes and spaces of a boundary of one region
classdef FEM_CouplingInterface < handle
    properties (Access = public)
        Region;                         % Region reference
        Boundaries;                     % Boundary numbers of region
        HarmonicsSin, HarmonicsCos;     % Used sin and cos values for coupling 
        HarmonicsAll;                   % Union of sin/cos for efficient rotation evaluation
        CouplingIndices;                % Indices for used sin/cos harmonics in all harmonics
        Radius;                         % Radius of Coupling Interface
        ArcLength;                      % Interface Length
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
        function obj = FEM_CouplingInterface(region, boundaryNr)
            obj.Region = region;
            obj.Boundaries = obj.Region.Boundaries(boundaryNr);
            obj.NumberBoundaryDOFs = numel(region.getBoundaryDOFs(boundaryNr));
            obj.generateSpacesAndMeshes();
        end

        function generateSpacesAndMeshes(obj)
            nodeRadius = vecnorm(obj.Region.Nodes(obj.Boundaries.dofs, :), 2, 2);
            obj.Radius = mean(nodeRadius);
            if any(abs(nodeRadius-obj.Radius) > 1e-5)
                warning('FEM_CouplingInterface:Radius','Coupling Interface Warning: Nodes might not lie on a circle! \t Min Radius: %.5f \t Max Radius: %.5f \n', min(nodeRadius), max(nodeRadius));
            end

            obj.ArcLength = 0;
            for iLine = 1:size(obj.Boundaries.lines, 1)
                % Length between two successive nodes
                L12 = vecnorm(obj.Region.Nodes(obj.Boundaries.lines(iLine,2),:) - obj.Region.Nodes(obj.Boundaries.lines(iLine,1),:));
                % Projection to circle
                obj.ArcLength = obj.ArcLength + 2*obj.Radius*asin(L12/(2*obj.Radius));
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
            nodes1 = obj.Boundaries.lines(:, 1);
            nodes2 = obj.Boundaries.lines(:, 2);
            t1 = atan2(obj.Region.Nodes(nodes1, 2), obj.Region.Nodes(nodes1, 1));
            t2 = atan2(obj.Region.Nodes(nodes2, 2), obj.Region.Nodes(nodes2, 1));
            theta1 = min(t1,t2);
            theta2 = max(t1,t2);
            % check for sign change at 180 degrees
            idxSgnChange = theta2-theta1>pi;
            if any(idxSgnChange)
                disp("Not yet confirmed to work for this, please confirm it works")
            end
            temp = theta1(idxSgnChange);
            theta1(idxSgnChange) = theta2(idxSgnChange);
            theta2(idxSgnChange) = temp+2*pi;

            % integration:
            n = obj.HarmonicsAll;
            N = numel(n);
            intN1sin_nt = (-n.*(theta1-theta2).*cos(n.*theta2) + sin(n.*theta1) - sin(n.*theta2))./(n.^2.*(theta1-theta2));
            intN2sin_nt = (n.*(theta1-theta2).*cos(n.*theta1) - sin(n.*theta1) + sin(n.*theta2))./(n.^2.*(theta1-theta2));
            intN1cos_nt = ((cos(n.*theta1) - cos(n.*theta2))./(theta1-theta2) + n.*sin(n.*theta2))./(n.^2);
            intN2cos_nt = (n.*(theta2-theta1).*sin(n.*theta1) - cos(n.*theta1) + cos(n.*theta2))./(n.^2.*(theta1-theta2));
            % correction for zero mode
            idxZero = find(obj.HarmonicsAll==0);
            if ~isempty(idxZero)
                intN1sin_nt(:, idxZero) = 0;
                intN2sin_nt(:, idxZero) = 0;
                intN1cos_nt(:, idxZero) = (theta2-theta1)/2;
                intN2cos_nt(:, idxZero) = (theta2-theta1)/2;
            end

            rows = [repmat(nodes1, N, 1); repmat(nodes2, N, 1); repmat(nodes1, N, 1); repmat(nodes2, N, 1)];
            nel = numel(nodes1);
            cols = [reshape(repmat(1:2:2*N, nel ,1), [], 1); ...
                    reshape(repmat(1:2:2*N, nel ,1), [], 1); ...
                    reshape(repmat(2:2:2*N, nel ,1), [], 1); ...
                    reshape(repmat(2:2:2*N, nel ,1), [], 1)];
            vals = [reshape(intN1sin_nt, [], 1); reshape(intN2sin_nt, [], 1); ...
                    reshape(intN1cos_nt, [], 1); reshape(intN2cos_nt, [], 1)]'*obj.Radius;
            obj.CouplingMatrixFull = sparse(rows, cols, vals, obj.Region.NumberDOF, 2*N);
            
            % OLD VERSION with loop
            % obj.CouplingMatrixFull = spalloc (obj.Region.NumberDOF, 2*numel(obj.HarmonicsAll), 2*numel(obj.HarmonicsAll)*numel(obj.Boundaries.dofs));
            % for iLine = 1:size(obj.Boundaries.lines, 1)
            %     node1 = obj.Boundaries.lines(iLine, 1);
            %     node2 = obj.Boundaries.lines(iLine, 2);
            %     t1 = atan2(obj.Region.Nodes(node1, 2), obj.Region.Nodes(node1, 1));
            %     t2 = atan2(obj.Region.Nodes(node2, 2), obj.Region.Nodes(node2, 1));
            %     theta1 = min(t1,t2);
            %     theta2 = max(t1,t2);
            %     % check for sign change at 180 degrees
            %     if theta2-theta1>pi
            %         temp = theta1;
            %         theta1 = theta2;
            %         theta2 = temp+2*pi;
            %     end
            %     % integration:
            %     n = obj.HarmonicsAll;
            %     intN1sin_nt = (-n.*(theta1-theta2).*cos(n*theta2) + sin(n*theta1) - sin(n*theta2))./(n.^2*(theta1-theta2));
            %     intN2sin_nt = (n.*(theta1-theta2).*cos(n*theta1) - sin(n*theta1) + sin(n*theta2))./(n.^2*(theta1-theta2));
            %     intN1cos_nt = ((cos(n*theta1) - cos(n*theta2))/(theta1-theta2) + n.*sin(n*theta2))./(n.^2);
            %     intN2cos_nt = (n.*(theta2-theta1).*sin(n*theta1) - cos(n*theta1) + cos(n*theta2))./(n.^2*(theta1-theta2));
            %     % correction for zero mode
            %     intN1sin_nt(obj.HarmonicsAll==0) = 0;
            %     intN2sin_nt(obj.HarmonicsAll==0) = 0;
            %     intN1cos_nt(obj.HarmonicsAll==0) = (theta2-theta1)/2;
            %     intN2cos_nt(obj.HarmonicsAll==0) = (theta2-theta1)/2;
            % 
            %     obj.CouplingMatrixFull(node1, 1:2:end) = obj.CouplingMatrixFull (node1, 1:2:end) + intN1sin_nt*obj.Radius;
            %     obj.CouplingMatrixFull(node2, 1:2:end) = obj.CouplingMatrixFull (node2, 1:2:end) + intN2sin_nt*obj.Radius;
            %     obj.CouplingMatrixFull(node1, 2:2:end) = obj.CouplingMatrixFull (node1, 2:2:end) + intN1cos_nt*obj.Radius;
            %     obj.CouplingMatrixFull(node2, 2:2:end) = obj.CouplingMatrixFull (node2, 2:2:end) + intN2cos_nt*obj.Radius;
            % end

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
            for iBnd = 1:numel(obj.Boundaries)
                xvals = NaN(3*size(obj.Boundaries(iBnd).lines, 1), 2);
                yvals = NaN(3*size(obj.Boundaries(iBnd).lines, 1), 2);
                for iLine = 1:size(obj.Boundaries(iBnd).lines, 1)
                    xvals((iLine-1)*3+1, 1) = obj.Region.NodesPlot(obj.Boundaries(iBnd).lines(iLine,1), 1);
                    xvals((iLine-1)*3+2, 1) = obj.Region.NodesPlot(obj.Boundaries(iBnd).lines(iLine,2), 1);
                    yvals((iLine-1)*3+1, 1) = obj.Region.NodesPlot(obj.Boundaries(iBnd).lines(iLine,1), 2);
                    yvals((iLine-1)*3+2, 1) = obj.Region.NodesPlot(obj.Boundaries(iBnd).lines(iLine,2), 2);
                end
                patch('XData',xvals,'YData',yvals,'EdgeColor', col, 'LineWidth', lw);
            end
        end
    end
end