classdef PLT_Lines < handle
    properties (Access = public)
        LinesPlot = struct();           % struct holding lines plot data and handles
        VectorPlot = struct();          % struct holding vector plot data and handles
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Line Plots %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeLines(obj, figureName, valueName, storeName, lineDimension, focusFigure)
            % create new figure if it does not exist
            figureObj = findobj('type', 'figure', 'name', figureName);
            figureName = matlab.lang.makeValidName(figureName);
            if ishandle(figureObj)
                obj.LinesPlot.(figureName).("Figure") = figureObj;
                if focusFigure
                    figure(figureObj); % focus figure
                end
            elseif ~isfield(obj.LinesPlot, figureName) || ~isfield(obj.LinesPlot.(figureName), "Figure") || ~ishandle(obj.LinesPlot.(figureName).("Figure"))
                obj.LinesPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
            end
            hold on
            % if line does not exist or plot dimension changed: reset values
            if ~isfield(obj.LinesPlot, figureName) || ~isfield(obj.LinesPlot.(figureName), valueName) || ~isfield(obj.LinesPlot.(figureName).(valueName), storeName)
                obj.LinesPlot.(figureName).(valueName).(storeName).Plot = gobjects(1);
                obj.resetLines(figureName, valueName, storeName, lineDimension)
            elseif any(obj.LinesPlot.(figureName).(valueName).(storeName).lineDimension ~= lineDimension, "all")
                obj.resetLines(figureName, valueName, storeName, lineDimension)
            end
        end

        function resetLines(obj, figureName, valueName, storeName, lineDimension)
            figureName = matlab.lang.makeValidName(figureName);
            obj.LinesPlot.(figureName).(valueName).(storeName).lineDimension = lineDimension;
            nLines = 0;
            for iDim = 1:size(lineDimension, 1)
                lCounter = lineDimension(iDim, 1)*3;
                obj.LinesPlot.(figureName).(valueName).(storeName).LineIndices{iDim, 1} = reshape((1:3:lCounter) + nLines, [], 1);
                obj.LinesPlot.(figureName).(valueName).(storeName).LineIndices{iDim, 2} = reshape((2:3:lCounter) + nLines, [], 1);
                nLines = nLines + lCounter;
            end
            obj.LinesPlot.(figureName).(valueName).(storeName).XData = NaN(nLines, 2);
            obj.LinesPlot.(figureName).(valueName).(storeName).YData = NaN(nLines, 2);
        end

        function setLines(obj, figureName, valueName, storeName, lineNumber, args)
            arguments
                obj;
                figureName="lines";
                valueName=figureName;
                storeName=valueName;
                lineNumber=1;
                args.XData (:,2) double = [];
                args.YData (:,2) double = [];
                args.ZData=[];
                args.CData=[];
            end
            figureName = matlab.lang.makeValidName(figureName);
            idxL1 = obj.LinesPlot.(figureName).(valueName).(storeName).LineIndices{lineNumber, 1};
            idxL2 = obj.LinesPlot.(figureName).(valueName).(storeName).LineIndices{lineNumber, 2};

            if ~isempty(args.XData)
                obj.LinesPlot.(figureName).(valueName).(storeName).XData(idxL1, 1) = reshape(args.XData(:,1), [], 1);
                obj.LinesPlot.(figureName).(valueName).(storeName).XData(idxL2, 1) = reshape(args.XData(:,2), [], 1);
            end
            if ~isempty(args.YData)
                obj.LinesPlot.(figureName).(valueName).(storeName).YData(idxL1, 1) = reshape(args.YData(:,1), [], 1);
                obj.LinesPlot.(figureName).(valueName).(storeName).YData(idxL2, 1) = reshape(args.YData(:,2), [], 1);
            end
        end

        function updateLines(obj, figureName, valueName, storeName, args)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                args.LineWidth = 2.5;
                args.EdgeColor = 'black';
                args.PickableParts {mustBeMember(args.PickableParts, ['none', 'visible', 'all'])} = 'none';
            end
            figureName = matlab.lang.makeValidName(figureName);
            if ~ishandle(obj.LinesPlot.(figureName).(valueName).(storeName).Plot)
                obj.LinesPlot.(figureName).(valueName).(storeName).Plot = patch('XData', obj.LinesPlot.(figureName).(valueName).(storeName).XData, ...
                    'YData', obj.LinesPlot.(figureName).(valueName).(storeName).YData, 'LineWidth', args.LineWidth, 'EdgeColor', args.EdgeColor, ...
                    'DeleteFcn', @(w, x) obj.removeDataLines(figureName, valueName, storeName), 'PickableParts', args.PickableParts);
            else
                set(obj.LinesPlot.(figureName).(valueName).(storeName).Plot, 'XData', obj.LinesPlot.(figureName).(valueName).(storeName).XData, ...
                    'YData', obj.LinesPlot.(figureName).(valueName).(storeName).YData, 'LineWidth', args.LineWidth, 'EdgeColor', args.EdgeColor);
            end
        end

        function removeDataLines(obj, figureName, valueName, storeName)
            figureName = matlab.lang.makeValidName(figureName);
            obj.LinesPlot.(figureName).(valueName) = rmfield(obj.LinesPlot.(figureName).(valueName), storeName);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Vector Plots %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeVector(obj, figureName, valueName, storeName, vectorDimension, focusFigure)
            % create new figure if it does not exist
            figureObj = findobj('type', 'figure', 'name', figureName);
            figureName = matlab.lang.makeValidName(figureName);
            if ishandle(figureObj)
                obj.VectorPlot.(figureName).("Figure") = figureObj;
                if focusFigure
                    figure(figureObj); % focus figure
                end
            elseif ~isfield(obj.VectorPlot, figureName) || ~isfield(obj.VectorPlot.(figureName), "Figure") || ~ishandle(obj.VectorPlot.(figureName).("Figure"))
                obj.VectorPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
            end
            hold on

            % if line does not exist or plot dimension changed: reset values
            if ~isfield(obj.VectorPlot, figureName) || ~isfield(obj.VectorPlot.(figureName), valueName) || ~isfield(obj.VectorPlot.(figureName).(valueName), storeName)
                obj.VectorPlot.(figureName).(valueName).(storeName).Plot = gobjects(1);
                obj.resetVector(figureName, valueName, storeName, vectorDimension)
            elseif any(obj.VectorPlot.(figureName).(valueName).(storeName).vectorDimension ~= vectorDimension, "all")
                obj.resetVector(figureName, valueName, storeName, vectorDimension)
            end
        end

        function resetVector(obj, figureName, valueName, storeName, vectorDimension)
            figureName = matlab.lang.makeValidName(figureName);
            obj.VectorPlot.(figureName).(valueName).(storeName).vectorDimension = vectorDimension;
            nVertices = 0;
            obj.VectorPlot.(figureName).(valueName).(storeName).VertexIndices = cell(size(vectorDimension, 1), 1);
            for iDim = 1:size(vectorDimension, 1)
                nVerticesPatch = vectorDimension(iDim, 1)*vectorDimension(iDim, 2);
                obj.VectorPlot.(figureName).(valueName).(storeName).VertexIndices{iDim, 1} = reshape((1:nVerticesPatch) + nVertices, [], 1);
                nVertices = nVertices + nVerticesPatch;
            end
            obj.VectorPlot.(figureName).(valueName).(storeName).XData = NaN(nVertices, 1);
            obj.VectorPlot.(figureName).(valueName).(storeName).YData = NaN(nVertices, 1);
            obj.VectorPlot.(figureName).(valueName).(storeName).UData = NaN(nVertices, 1);
            obj.VectorPlot.(figureName).(valueName).(storeName).VData = NaN(nVertices, 1);
        end

        function setVector(obj, figureName, valueName, storeName, patchNumber, args)
            arguments
                obj;
                figureName="vector";
                valueName=figureName;
                storeName=valueName;
                patchNumber=1;
                args.XData=[];
                args.YData=[];
                args.UData=[];
                args.VData=[];
                args.Rotation=0;
            end
            figureName = matlab.lang.makeValidName(figureName);
            idx = obj.VectorPlot.(figureName).(valueName).(storeName).VertexIndices{patchNumber};
            if ~isempty(args.XData)
                obj.VectorPlot.(figureName).(valueName).(storeName).XData(idx) = reshape(args.XData, [], 1);
            end
            if ~isempty(args.YData)
                obj.VectorPlot.(figureName).(valueName).(storeName).YData(idx) = reshape(args.YData, [], 1);
            end
            if args.Rotation~=0
                assert(~isempty(args.UData) && ~isempty(args.VData));
                R = [cos(args.Rotation), sin(args.Rotation); -sin(args.Rotation), cos(args.Rotation)];
                rotData = [reshape(args.UData, [], 1), reshape(args.VData, [],1)] *R;
                args.UData = rotData(:,1);
                args.VData = rotData(:,2);
            end
            if ~isempty(args.UData)
                obj.VectorPlot.(figureName).(valueName).(storeName).UData(idx) = reshape(args.UData, [], 1);
            end
            if ~isempty(args.VData)
                obj.VectorPlot.(figureName).(valueName).(storeName).VData(idx) = reshape(args.VData, [], 1);
            end
        end

        function updateVector(obj, figureName, valueName, storeName, args)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                args.Color = "black";
                args.Scaling (1,1) double = 0.5;    % only set initially and not updated
                args.LineWidth (1,1) double = 0.5;
            end
            figureName = matlab.lang.makeValidName(figureName);
            
            x = obj.VectorPlot.(figureName).(valueName).(storeName).XData;
            y = obj.VectorPlot.(figureName).(valueName).(storeName).YData;
            u = obj.VectorPlot.(figureName).(valueName).(storeName).UData;
            v = obj.VectorPlot.(figureName).(valueName).(storeName).VData;
            ind = find((u.^2+v.^2).^0.5);
            
            if ~ishandle(obj.VectorPlot.(figureName).(valueName).(storeName).Plot)
                obj.VectorPlot.(figureName).(valueName).(storeName).Plot = ...
                quiver(x(ind), y(ind), u(ind), v(ind), ...
                       args.Scaling, 'color', args.Color, 'LineWidth', args.LineWidth, ...
                       'DeleteFcn', @(w, x) obj.removeDataVector(figureName, valueName, storeName));
                % uistack(obj.VectorPlot.(figureName).(valueName).(storeName).Plot, 'top');
            else
                set(obj.VectorPlot.(figureName).(valueName).(storeName).Plot, ...
                    'XData', x(ind), 'YData', y(ind), 'UData', u(ind), 'VData', v(ind), ...
                    'Color', args.Color, 'LineWidth', args.LineWidth);
            end
        end

        function removeDataVector(obj, figureName, valueName, storeName)
            figureName = matlab.lang.makeValidName(figureName);
            obj.VectorPlot.(figureName).(valueName) = rmfield(obj.VectorPlot.(figureName).(valueName), storeName);
        end

    end
end
