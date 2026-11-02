classdef PLT_TriangularPattern < PLT_Lines
    properties (Access = public)
        ScalarPlot = struct();          % struct holding scalar plot data and handles
        ColorPlot = struct();           % struct holding color plot data and handles
        ContourPlot = struct();         % struct holding contour plot data and handles
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Scalar Plots %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeScalar(obj, figureName, valueName, storeName, triangleDimensions, focusFigure)
            % create new figure if it does not exist
            figureObj = findobj('type', 'figure', 'name', figureName);
            if ishandle(figureObj)
                obj.ScalarPlot.(figureName).("Figure") = figureObj;
                if focusFigure
                    figure(figureObj); % focus figure
                end
            elseif ~isfield(obj.ScalarPlot, figureName) || ~isfield(obj.ScalarPlot.(figureName), "Figure") || ~ishandle(obj.ScalarPlot.(figureName).("Figure"))
                obj.ScalarPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
            end
            hold on
            if ~isfield(obj.ScalarPlot, figureName) || ~isfield(obj.ScalarPlot.(figureName), valueName) || ~isfield(obj.ScalarPlot.(figureName).(valueName), storeName)
                obj.ScalarPlot.(figureName).(valueName).(storeName).Plot = gobjects(1);
                obj.resetScalar(figureName, valueName, storeName, triangleDimensions);
            elseif any(obj.ScalarPlot.(figureName).(valueName).(storeName).triangleDimension ~= triangleDimensions, "all")
                obj.resetScalar(figureName, valueName, storeName, triangleDimensions);
            end
        end

        function resetScalar(obj, figureName, valueName, storeName, triangleDimensions)
            obj.ScalarPlot.(figureName).(valueName).(storeName).triangleDimension = triangleDimensions;
            nFaces = 0;
            nVertices = 0;
            for iDim = 1:size(triangleDimensions, 1)
                nVerticesTriangle = triangleDimensions(iDim, 1)*triangleDimensions(iDim, 2);
                obj.ScalarPlot.(figureName).(valueName).(storeName).VertexIndices{iDim, 1} = reshape((1:nVerticesTriangle) + nVertices, [], 1);
                nVertices = nVertices + nVerticesTriangle;
                nFaces = nFaces + triangleDimensions(iDim, 1);
            end
            obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices = NaN(nVertices, 3);
            obj.ScalarPlot.(figureName).(valueName).(storeName).Colors = NaN(nVertices, 1);
            obj.ScalarPlot.(figureName).(valueName).(storeName).Faces = NaN(nFaces, 3);
            % directly update face indices
            nFaces = 0;
            nVertices = 0;
            for iDim = 1:size(triangleDimensions, 1)
                nFacesTriangle = triangleDimensions(iDim, 1);
                nVerticesTriangle = triangleDimensions(iDim, 1)*triangleDimensions(iDim, 2);
                % asdf = triangleDimensions(iDim, 1);
                obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 1) = (1:3:nVerticesTriangle) + nVertices;
                obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 2) = (2:3:nVerticesTriangle) + nVertices;
                obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 3) = (3:3:nVerticesTriangle) + nVertices;
                % obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 1) = (1:asdf) + nVertices;
                % obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 2) = (asdf+1:2*asdf) + nVertices;
                % obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:nFacesTriangle) + nFaces, 3) = (2*asdf+1:3*asdf) + nVertices;
                nFaces = nFaces + nFacesTriangle;
                nVertices = nVertices + nVerticesTriangle;
            end
        end

        function setScalar(obj, figureName, valueName, storeName, patchNumber, args)
            arguments
                obj;
                figureName="scalar";
                valueName=figureName;
                storeName=valueName;
                patchNumber=1;
                args.XData=[];
                args.YData=[];
                args.ZData=[];
                args.CData=[];
                args.flat (1,1) logical = false;
                args.absolute (1,1) logical = false;
                args.scaling = 1;
            end
            idx = obj.ScalarPlot.(figureName).(valueName).(storeName).VertexIndices{patchNumber};
            if ~isempty(args.XData)
                obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices(idx, 1) = reshape(args.XData, [], 1);
            end
            if ~isempty(args.YData)
                obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices(idx, 2) = reshape(args.YData, [], 1);
            end
            if ~isempty(args.ZData)
                if args.flat
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices(idx, 3) = zeros(size(reshape(args.ZData, [], 1)));
                elseif args.absolute
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices(idx, 3) = args.scaling*reshape(abs(args.ZData), [], 1);
                else
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices(idx, 3) = args.scaling*reshape(args.ZData, [], 1);
                end
            end
            if ~isempty(args.CData)
                if args.absolute
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Colors(idx) = args.scaling*reshape(abs(args.CData), [], 1);
                else
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Colors(idx) = args.scaling*reshape(args.CData, [], 1);
                end
            end
        end

        function updateScalar(obj, figureName, valueName, storeName, args)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                args.Shading = 'interp';
                args.EdgeColor = 'none';
            end
            if ~ishandle(obj.ScalarPlot.(figureName).(valueName).(storeName).Plot)
                obj.ScalarPlot.(figureName).(valueName).(storeName).Plot = patch('Faces', obj.ScalarPlot.(figureName).(valueName).(storeName).Faces, 'Vertices', obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices,...
                    'FaceVertexCData', obj.ScalarPlot.(figureName).(valueName).(storeName).Colors, 'FaceColor', args.Shading, 'EdgeColor', args.EdgeColor, ...
                    'DeleteFcn', @(w, x) obj.removeDataScalar(figureName, valueName, storeName));
                % necessary for initialization of the datatip template
                datatip(obj.ScalarPlot.(figureName).(valueName).(storeName).Plot, 'Visible', 'off');
                obj.ScalarPlot.(figureName).(valueName).(storeName).Plot.DataTipTemplate.DataTipRows(3,1).Value = 'CData';
            else
                set(obj.ScalarPlot.(figureName).(valueName).(storeName).Plot, 'Faces', obj.ScalarPlot.(figureName).(valueName).(storeName).Faces, 'Vertices', obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices,...
                    'FaceVertexCData', obj.ScalarPlot.(figureName).(valueName).(storeName).Colors, 'FaceColor', args.Shading, 'EdgeColor', args.EdgeColor);
            end
        end

        function removeDataScalar(obj, figureName, valueName, storeName)
            obj.ScalarPlot.(figureName).(valueName) = rmfield(obj.ScalarPlot.(figureName).(valueName), storeName);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Color Plots %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeColor(obj, figureName, valueName, storeName, numberRegions, focusFigure)
            % create new figure if it does not exist
            figureObj = findobj('type', 'figure', 'name', figureName);
            figureName = matlab.lang.makeValidName(figureName);
            if ishandle(figureObj)
                obj.ColorPlot.(figureName).("Figure") = figureObj;
                if focusFigure
                    figure(figureObj); % focus figure
                end
            elseif ~isfield(obj.ColorPlot, figureName) || ~isfield(obj.ColorPlot.(figureName), "Figure") || ~ishandle(obj.ColorPlot.(figureName).("Figure"))
                obj.ColorPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
            end
            hold on
            % if field does not exist or plot dimension changed: reset values
            if ~isfield(obj.ColorPlot, figureName) || ~isfield(obj.ColorPlot.(figureName), valueName) || ~isfield(obj.ColorPlot.(figureName).(valueName), storeName)
                obj.ColorPlot.(figureName).(valueName).(storeName).Plot = gobjects(numberRegions, 1);
                obj.resetColor(figureName, valueName, storeName, numberRegions);
            elseif any(obj.ColorPlot.(figureName).(valueName).(storeName).NumberRegions ~= numberRegions, "all")
                obj.resetColor(figureName, valueName, storeName, numberRegions);
            end
        end

        function resetColor(obj, figureName, valueName, storeName, numberRegions)
            obj.ColorPlot.(figureName).(valueName).(storeName).NumberRegions = numberRegions;
            % for iDim = 1:size(numberRegions, 1)
            %     obj.ColorPlot.(figureName).(valueName).(storeName).Faces = NaN(numberRegions(iDim,:));
            %     obj.ColorPlot.(figureName).(valueName).(storeName).Vertices = NaN(numberRegions(iDim,:));
            % end
            obj.ColorPlot.(figureName).(valueName).(storeName).Color = cell(numberRegions, 1);
            obj.ColorPlot.(figureName).(valueName).(storeName).Alpha = cell(numberRegions, 1);
        end

        function setColor(obj, figureName, valueName, storeName, regionNumber, args)
            arguments
                obj;
                figureName="color";
                valueName=figureName;
                storeName=valueName;
                regionNumber=1;
                args.Faces=[];
                args.Vertices=[];
                args.Color="white";
                args.Alpha=1;
            end

            if ~isempty(args.Faces)
                obj.ColorPlot.(figureName).(valueName).(storeName).Faces{regionNumber} = args.Faces;
            end
            if ~isempty(args.Vertices)
                obj.ColorPlot.(figureName).(valueName).(storeName).Vertices{regionNumber} = args.Vertices;
            end
            obj.ColorPlot.(figureName).(valueName).(storeName).Color{regionNumber} = args.Color;
            obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{regionNumber} = args.Alpha;
        end

        function updateColor(obj, figureName, valueName, storeName, args)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                args.EdgeColor = "none";
                args.LineWidth = 1.0;
                % drawOrder: when exporting to pdf, the order matters, as annoying lines appear when switching the color between patches
                args.drawOrder = 1:obj.ColorPlot.(figureName).(valueName).(storeName).NumberRegions;
            end

            for region = args.drawOrder
                if ~ishandle(obj.ColorPlot.(figureName).(valueName).(storeName).Plot(region))
                    obj.ColorPlot.(figureName).(valueName).(storeName).Plot(region) = ...
                        trisurf(obj.ColorPlot.(figureName).(valueName).(storeName).Faces{region}, ...
                        obj.ColorPlot.(figureName).(valueName).(storeName).Vertices{region}(:,1), ...
                        obj.ColorPlot.(figureName).(valueName).(storeName).Vertices{region}(:,2), ...
                        zeros(size(obj.ColorPlot.(figureName).(valueName).(storeName).Vertices{region}(:,1))), ...
                        'FaceColor', obj.ColorPlot.(figureName).(valueName).(storeName).Color{region}, ...
                        'FaceAlpha', obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{region}, ...
                        'EdgeColor', args.EdgeColor, 'LineWidth', args.LineWidth, ...
                        'DeleteFcn', @(w, x) obj.removeDataColor(figureName, valueName, storeName));
                else
                    set(obj.ColorPlot.(figureName).(valueName).(storeName).Plot(region), ...
                        'Faces', obj.ColorPlot.(figureName).(valueName).(storeName).Faces{region}, ...
                        'Vertices', obj.ColorPlot.(figureName).(valueName).(storeName).Vertices{region},...
                        'FaceColor', obj.ColorPlot.(figureName).(valueName).(storeName).Color{region}, ...
                        'FaceAlpha', obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{region}, ...
                        'EdgeColor', args.EdgeColor, 'LineWidth', args.LineWidth);
                end
            end
        end

        function removeDataColor(obj, figureName, valueName, storeName)
            if isfield(obj.ColorPlot.(figureName).(valueName), storeName)
                obj.ColorPlot.(figureName).(valueName) = rmfield(obj.ColorPlot.(figureName).(valueName), storeName);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Contour Plots %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function initializeContour(obj, figureName, valueName, storeName, patchDimensions)
        %     % create new figure if it does not exist
        %     figureObj = findobj('type', 'figure', 'name', figureName);
        %     if ishandle(figureObj)
        %         obj.ContourPlot.(figureName).("Figure") = figureObj;
        %         figure(figureObj);
        %     elseif ~isfield(obj.ContourPlot, figureName) || ~isfield(obj.ContourPlot.(figureName), "Figure") || ~ishandle(obj.ContourPlot.(figureName).("Figure"))
        %         obj.ContourPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
        %     end
        %     hold on
        %     obj.ContourPlot.(figureName).(valueName).(storeName).Plot = gobjects(size(patchDimensions, 1), 1);
        %     obj.ContourPlot.(figureName).(valueName).(storeName).patchDimension = patchDimensions;
        % 
        %     for iDim = 1:size(patchDimensions, 1)
        %         obj.ContourPlot.(figureName).(valueName).(storeName).XData{iDim, 1} = NaN(patchDimensions(iDim,:));
        %         obj.ContourPlot.(figureName).(valueName).(storeName).YData{iDim, 1} = NaN(patchDimensions(iDim,:));
        %         obj.ContourPlot.(figureName).(valueName).(storeName).ZData{iDim, 1} = NaN(patchDimensions(iDim,:));
        %     end
        % end
        % 
        % function setContour(obj, figureName, valueName, storeName, patchNumber, args)
        %     arguments
        %         obj;
        %         figureName="contour";
        %         valueName=figureName;
        %         storeName=valueName;
        %         patchNumber=1;
        %         args.XData=[];
        %         args.YData=[];
        %         args.ZData=[];
        %         args.absolute (1,1) logical = false;
        %         args.scaling = 1;
        %     end
        % 
        %     if ~isempty(args.XData)
        %         obj.ContourPlot.(figureName).(valueName).(storeName).XData{patchNumber} = args.XData;
        %     end
        %     if ~isempty(args.YData)
        %         obj.ContourPlot.(figureName).(valueName).(storeName).YData{patchNumber} = args.YData;
        %     end
        %     if ~isempty(args.ZData)
        %         if args.absolute
        %             obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patchNumber} = args.scaling*abs(args.ZData);
        %         else
        %             obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patchNumber} = args.scaling*args.ZData;
        %         end
        %     end
        % end
        % 
        % function updateContour(obj, figureName, valueName, storeName, args)
        %     arguments
        %         obj;
        %         figureName;
        %         valueName;
        %         storeName;
        %         args.LevelLimits = [];
        %         args.nLevels (1,1) double = 20;
        %         args.LineColor = "black";
        %         args.LineWidth (1,1) double = 0.5;
        %         args.filled (1,1) logical = true;
        %         % drawOrder: when exporting to pdf, the order matters, as annoying lines appear when switching the color between patches
        %         args.drawOrder = 1:size(obj.ContourPlot.(figureName).(valueName).(storeName).patchDimension, 1);
        %     end
        %     if isempty(args.LevelLimits)
        %         minVal = min(cellfun(@(x) min(x, [], "all"), obj.ContourPlot.(figureName).(valueName).(storeName).ZData));
        %         maxVal = max(cellfun(@(x) max(x, [], "all"), obj.ContourPlot.(figureName).(valueName).(storeName).ZData));
        %         if abs(maxVal-minVal) < 1e-10
        %             maxVal = maxVal + 1e-10;
        %         end
        %         levels = linspace(minVal, maxVal, args.nLevels);
        %     else
        %        levels = linspace(args.LevelLimits(1), args.LevelLimits(2), args.nLevels);
        %     end
        % 
        %     for patch = args.drawOrder
        %         if ~ishandle(obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch))
        %             if args.filled
        %                 [~, obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch)] = ...
        %                     contourf(obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
        %                              obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
        %                              obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch}, ...
        %                              'LevelList', levels, 'LineWidth', args.LineWidth, ...
        %                              'LineColor', args.LineColor, 'DeleteFcn', @(w, x) obj.removeDataContour(figureName, valueName, storeName));
        % 
        %             else
        %                 [~, obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch)] = ...
        %                     contour(obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
        %                             obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
        %                             obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch},...
        %                             'LevelList', levels, 'LineWidth', args.LineWidth, ...
        %                             'LineColor', args.LineColor, 'DeleteFcn', @(w, x) obj.removeDataContour(figureName, valueName, storeName));
        %             end
        %         else
        %             set(obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch), ...
        %                 'XData', obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
        %                 'YData', obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
        %                 'ZData', obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch}, ...
        %                 'LevelList', levels, 'LineWidth', args.LineWidth, ...
        %                 'LineColor', args.LineColor);
        %         end
        %     end
        %     if args.filled
        %         clim([levels(1), levels(end)]);
        %     end
        % end
        % 
        % function removeDataContour(obj, figureName, valueName, storeName)
        %     if isfield(obj.ContourPlot.(figureName).(valueName), storeName)
        %         obj.ContourPlot.(figureName).(valueName) = rmfield(obj.ContourPlot.(figureName).(valueName), storeName);
        %     end
        % end

    end
end
