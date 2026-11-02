classdef PLT_RectangularPattern < PLT_Lines
    properties (Access = public)
        ScalarPlot = struct();          % struct holding scalar plot data and handles
        ColorPlot = struct();           % struct holding color plot data and handles
        ContourPlot = struct();         % struct holding contour plot data and handles
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Scalar Plots %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeScalar(obj, figureName, valueName, storeName, patchDimensions, focusFigure)
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
                obj.resetScalar(figureName, valueName, storeName, patchDimensions)
            elseif any(obj.ScalarPlot.(figureName).(valueName).(storeName).patchDimension ~= patchDimensions, "all")
                obj.resetScalar(figureName, valueName, storeName, patchDimensions)
            end
        end

        function resetScalar(obj, figureName, valueName, storeName, patchDimensions)
            obj.ScalarPlot.(figureName).(valueName).(storeName).patchDimension = patchDimensions;
            nFaces = 0;
            nVertices = 0;
            for iDim = 1:size(patchDimensions, 1)
                nVerticesPatch = patchDimensions(iDim, 1)*patchDimensions(iDim, 2);
                obj.ScalarPlot.(figureName).(valueName).(storeName).VertexIndices{iDim, 1} = reshape((1:nVerticesPatch) + nVertices, [], 1);
                nVertices = nVertices + nVerticesPatch;
                nFaces = nFaces + (patchDimensions(iDim, 1)-1)*(patchDimensions(iDim, 2)-1);
            end
            obj.ScalarPlot.(figureName).(valueName).(storeName).Vertices = NaN(nVertices, 3);
            obj.ScalarPlot.(figureName).(valueName).(storeName).Colors = NaN(nVertices, 1);
            obj.ScalarPlot.(figureName).(valueName).(storeName).Faces = NaN(nFaces, 4);
            % directly update face indices
            nFaces = 0;
            nVertices = 0;
            for iDim = 1:size(patchDimensions, 1)
                nx = patchDimensions(iDim, 1);
                ny = patchDimensions(iDim, 2);
                for iny = 0:ny-2
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:(nx-1)) + iny*(nx-1) + nFaces, 1) = (1:(nx-1)) + nx*iny     + nVertices;
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:(nx-1)) + iny*(nx-1) + nFaces, 2) = (2:nx)     + nx*iny     + nVertices;
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:(nx-1)) + iny*(nx-1) + nFaces, 3) = (2:nx)     + nx*(iny+1) + nVertices;
                    obj.ScalarPlot.(figureName).(valueName).(storeName).Faces((1:(nx-1)) + iny*(nx-1) + nFaces, 4) = (1:(nx-1)) + nx*(iny+1) + nVertices;
                end
                nFacesPatch = (patchDimensions(iDim, 1)-1)*(patchDimensions(iDim, 2)-1);
                nVerticesPatch = patchDimensions(iDim, 1)*patchDimensions(iDim, 2);
                nFaces = nFaces + nFacesPatch;
                nVertices = nVertices + nVerticesPatch;
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
        function initializeColor(obj, figureName, valueName, storeName, patchDimensions, focusFigure)
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
                obj.ColorPlot.(figureName).(valueName).(storeName).Plot = gobjects(size(patchDimensions, 1), 1);
                obj.resetColor(figureName, valueName, storeName, patchDimensions);
            elseif any(obj.ColorPlot.(figureName).(valueName).(storeName).patchDimension ~= patchDimensions, "all")
                obj.resetColor(figureName, valueName, storeName, patchDimensions);
            end
        end

        function resetColor(obj, figureName, valueName, storeName, patchDimensions)
            figureName = matlab.lang.makeValidName(figureName);
            obj.ColorPlot.(figureName).(valueName).(storeName).patchDimension = patchDimensions;
            for iDim = 1:size(patchDimensions, 1)
                obj.ColorPlot.(figureName).(valueName).(storeName).XData{iDim, 1} = NaN(patchDimensions(iDim,:));
                obj.ColorPlot.(figureName).(valueName).(storeName).YData{iDim, 1} = NaN(patchDimensions(iDim,:));
            end
            obj.ColorPlot.(figureName).(valueName).(storeName).Color = cell(size(patchDimensions, 1), 1);
            obj.ColorPlot.(figureName).(valueName).(storeName).Alpha = cell(size(patchDimensions, 1), 1);
        end

        function setColor(obj, figureName, valueName, storeName, patchNumber, args)
            arguments
                obj;
                figureName="color";
                valueName=figureName;
                storeName=valueName;
                patchNumber=1;
                args.XData=[];
                args.YData=[];
                args.Color="white";
                args.Alpha=1;
            end
            figureName = matlab.lang.makeValidName(figureName);
            if ~isempty(args.XData)
                obj.ColorPlot.(figureName).(valueName).(storeName).XData{patchNumber} = args.XData;
            end
            if ~isempty(args.YData)
                obj.ColorPlot.(figureName).(valueName).(storeName).YData{patchNumber} = args.YData;
            end
            obj.ColorPlot.(figureName).(valueName).(storeName).Color{patchNumber} = args.Color;
            obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{patchNumber} = args.Alpha;
        end

        function updateColor(obj, figureName, valueName, storeName, drawOrder)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                % drawOrder: when exporting to pdf, the order matters, as annoying lines appear when switching the color between patches
                drawOrder = 1:size(obj.ColorPlot.(figureName).(valueName).(storeName).patchDimension, 1);
            end
            figureName = matlab.lang.makeValidName(figureName);
            for patch = drawOrder
                if ~ishandle(obj.ColorPlot.(figureName).(valueName).(storeName).Plot(patch))
                    obj.ColorPlot.(figureName).(valueName).(storeName).Plot(patch) = ...
                        surf(obj.ColorPlot.(figureName).(valueName).(storeName).XData{patch}, ...
                        obj.ColorPlot.(figureName).(valueName).(storeName).YData{patch}, ...
                        zeros(size(obj.ColorPlot.(figureName).(valueName).(storeName).XData{patch})), ...
                        'FaceColor', obj.ColorPlot.(figureName).(valueName).(storeName).Color{patch}, ...
                        'FaceAlpha', obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{patch}, ...
                        'EdgeColor', 'none', 'DeleteFcn', @(w, x) obj.removeDataColor(figureName, valueName, storeName));
                else
                    set(obj.ColorPlot.(figureName).(valueName).(storeName).Plot(patch), ...
                        'XData', obj.ColorPlot.(figureName).(valueName).(storeName).XData{patch}, ...
                        'YData', obj.ColorPlot.(figureName).(valueName).(storeName).YData{patch},...
                        'ZData', zeros(size(obj.ColorPlot.(figureName).(valueName).(storeName).XData{patch})), ...
                        'FaceColor', obj.ColorPlot.(figureName).(valueName).(storeName).Color{patch}, ...
                        'FaceAlpha', obj.ColorPlot.(figureName).(valueName).(storeName).Alpha{patch});
                end
            end
        end

        function removeDataColor(obj, figureName, valueName, storeName)
            figureName = matlab.lang.makeValidName(figureName);
            if isfield(obj.ColorPlot.(figureName).(valueName), storeName)
                obj.ColorPlot.(figureName).(valueName) = rmfield(obj.ColorPlot.(figureName).(valueName), storeName);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Contour Plots %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeContour(obj, figureName, valueName, storeName, patchDimensions, focusFigure)
            % create new figure if it does not exist
            figureObj = findobj('type', 'figure', 'name', figureName);
            if ishandle(figureObj)
                obj.ContourPlot.(figureName).("Figure") = figureObj;
                if focusFigure
                    figure(figureObj); % focus figure
                end
            elseif ~isfield(obj.ContourPlot, figureName) || ~isfield(obj.ContourPlot.(figureName), "Figure") || ~ishandle(obj.ContourPlot.(figureName).("Figure"))
                obj.ContourPlot.(figureName).("Figure") = figure('Name', figureName, 'NumberTitle', 'off');
            end
            hold on

            hold on
            if ~isfield(obj.ContourPlot, figureName) || ~isfield(obj.ContourPlot.(figureName), valueName) || ~isfield(obj.ContourPlot.(figureName).(valueName), storeName)
                obj.ContourPlot.(figureName).(valueName).(storeName).Plot = gobjects(size(patchDimensions, 1), 1);
                obj.resetContour(figureName, valueName, storeName, patchDimensions)
            elseif any(obj.ContourPlot.(figureName).(valueName).(storeName).patchDimension ~= patchDimensions, "all")
                obj.resetContour(figureName, valueName, storeName, patchDimensions)
            end
        end

        function resetContour(obj, figureName, valueName, storeName, patchDimensions)
            obj.ContourPlot.(figureName).(valueName).(storeName).patchDimension = patchDimensions;

            for iDim = 1:size(patchDimensions, 1)
                obj.ContourPlot.(figureName).(valueName).(storeName).XData{iDim, 1} = NaN(patchDimensions(iDim,:));
                obj.ContourPlot.(figureName).(valueName).(storeName).YData{iDim, 1} = NaN(patchDimensions(iDim,:));
                obj.ContourPlot.(figureName).(valueName).(storeName).ZData{iDim, 1} = NaN(patchDimensions(iDim,:));
            end
        end

        function setContour(obj, figureName, valueName, storeName, patchNumber, args)
            arguments
                obj;
                figureName="contour";
                valueName=figureName;
                storeName=valueName;
                patchNumber=1;
                args.XData=[];
                args.YData=[];
                args.ZData=[];
                args.absolute (1,1) logical = false;
                args.scaling = 1;
            end

            if ~isempty(args.XData)
                obj.ContourPlot.(figureName).(valueName).(storeName).XData{patchNumber} = args.XData;
            end
            if ~isempty(args.YData)
                obj.ContourPlot.(figureName).(valueName).(storeName).YData{patchNumber} = args.YData;
            end
            if ~isempty(args.ZData)
                if args.absolute
                    obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patchNumber} = args.scaling*abs(args.ZData);
                else
                    obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patchNumber} = args.scaling*args.ZData;
                end
            end
        end

        function updateContour(obj, figureName, valueName, storeName, args)
            arguments
                obj;
                figureName;
                valueName;
                storeName;
                args.LevelLimits = [];
                args.nLevels (1,1) double = 20;
                args.LineColor = "black";
                args.LineWidth (1,1) double = 0.5;
                args.filled (1,1) logical = true;
                % drawOrder: when exporting to pdf, the order matters, as annoying lines appear when switching the color between patches
                args.drawOrder = 1:size(obj.ContourPlot.(figureName).(valueName).(storeName).patchDimension, 1);
            end
            if isempty(args.LevelLimits)
                minVal = min(cellfun(@(x) min(x, [], "all"), obj.ContourPlot.(figureName).(valueName).(storeName).ZData));
                maxVal = max(cellfun(@(x) max(x, [], "all"), obj.ContourPlot.(figureName).(valueName).(storeName).ZData));
                if abs(maxVal-minVal) < 1e-10
                    maxVal = maxVal + 1e-10;
                end
                levels = linspace(minVal, maxVal, args.nLevels);
            else
               levels = linspace(args.LevelLimits(1), args.LevelLimits(2), args.nLevels);
            end

            for patch = args.drawOrder
                if ~ishandle(obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch))
                    if args.filled
                        [~, obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch)] = ...
                            contourf(obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
                                     obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
                                     obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch}, ...
                                     'LevelList', levels, 'LineWidth', args.LineWidth, ...
                                     'LineColor', args.LineColor, 'DeleteFcn', @(w, x) obj.removeDataContour(figureName, valueName, storeName));

                    else
                        [~, obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch)] = ...
                            contour(obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
                                    obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
                                    obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch},...
                                    'LevelList', levels, 'LineWidth', args.LineWidth, ...
                                    'LineColor', args.LineColor, 'DeleteFcn', @(w, x) obj.removeDataContour(figureName, valueName, storeName));
                    end
                else
                    set(obj.ContourPlot.(figureName).(valueName).(storeName).Plot(patch), ...
                        'XData', obj.ContourPlot.(figureName).(valueName).(storeName).XData{patch}, ...
                        'YData', obj.ContourPlot.(figureName).(valueName).(storeName).YData{patch}, ...
                        'ZData', obj.ContourPlot.(figureName).(valueName).(storeName).ZData{patch}, ...
                        'LevelList', levels, 'LineWidth', args.LineWidth, ...
                        'LineColor', args.LineColor);
                end
            end
            if args.filled
                clim([levels(1), levels(end)]);
            end
        end

        function removeDataContour(obj, figureName, valueName, storeName)
            if isfield(obj.ContourPlot.(figureName).(valueName), storeName)
                obj.ContourPlot.(figureName).(valueName) = rmfield(obj.ContourPlot.(figureName).(valueName), storeName);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Exports Plots %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function exportParaview(obj, figureName, time)
            %  TBD go through struct with name figureName, see which types it is
            % and export one after another to files
        end

        % function exportParaviewAz(obj, filename, valueIndex, valueList)
        %     if ~isfolder("paraview")
        %         mkdir("paraview")
        %     end
        %     if ~exist("filename", "var")
        %         filename = 'paraview/Az';
        %     end
        %     if nargin <= 2
        %         valueIndex = 1;
        %         valueList = 1;
        %     end
        %     obj.checkPlotValues();
        %     obj.checkRotationUpdate();
        % 
        %     str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
        %     str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
        %     str3 = cat (2, '</Collection>\n', '</VTKFile> \n');
        % 
        %     if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
        %         pvd_filename = cat (2, filename, '.pvd');
        %     else
        %         pvd_filename = filename;
        %         filename = filename (1:end-4);
        %     end
        % 
        %     if valueIndex == 1
        %         fid = fopen (pvd_filename, 'w');
        %         fprintf (fid, str1);
        %     else
        %         fid = fopen (pvd_filename, 'a');
        %     end
        % 
        %     if (fid < 0)
        %         error ('exportSolutionParaView: could not open file %s', pvd_filename);
        %     end
        % 
        %     ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
        %     if (isempty (ind)); ind = 1; end
        % 
        %     for iPatch = 1:obj.NumberPatches
        %         filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(valueIndex-1), '_', num2str (iPatch));
        %         filename_patch = cat (2, filename, '_', num2str(valueIndex-1), '_', num2str (iPatch));
        %         vts_pts = [];
        % 
        %         vts_pts(1, :, :) = obj.Geometry(iPatch).PlotPoints{1};
        %         vts_pts(2, :, :) = obj.Geometry(iPatch).PlotPoints{2};
        %         Az_patch = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
        %         Az_points = sum(obj.Geometry(iPatch).PlotValues{1}.* reshape(Az_patch, 1, 1, []), 3) / obj.Length;
        % 
        %         msh_to_vtk (vts_pts, Az_points, filename_patch, 'A_z');
        % 
        %         fprintf (fid, str2, valueList(valueIndex), iPatch, filename_patch_without_path);
        %     end
        % 
        %     if valueIndex == numel(valueList)
        %         fprintf (fid, str3);
        %     end
        %     fclose (fid);
        % end
        % 
        % function exportTransientParaview(obj, filename, exportTime)
        %     if ~exist("exportTime", "var")
        %         exportTime = obj.Solution.Time;
        %     end
        %     progressbar('Exporting to Paraview: ')
        %     for ti = 1:numel(exportTime)
        %         obj.MagneticPotential = interp1(obj.Solution.Time, obj.Solution.Values', exportTime(ti))';
        %         obj.exportParaviewAz(['paraview/' filename '_Az'], ti, exportTime);
        %         obj.exportParaviewB(['paraview/' filename '_B'], ti, exportTime);
        %         progressbar(ti/numel(exportTime));
        %     end
        %     progressbar('finished');
        % end
        % 
        % % Export Magnetic B Field to ParaView
        % function exportParaviewB(obj, filename, currentIndex, dataList)
        %     if ~isfolder("paraview")
        %         mkdir("paraview")
        %     end
        %     if ~exist("filename", "var")
        %         filename = 'paraview/B';
        %     end
        %     if nargin <= 2
        %         currentIndex = 1;
        %         dataList = 1;
        %     end
        %     obj.checkPlotValues();
        % 
        %     str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
        %     str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
        %     str3 = cat (2, '</Collection>\n', '</VTKFile> \n');
        % 
        %     if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
        %         pvd_filename = cat (2, filename, '.pvd');
        %     else
        %         pvd_filename = filename;
        %         filename = filename (1:end-4);
        %     end
        % 
        %     if currentIndex == 1
        %         fid = fopen (pvd_filename, 'w');
        %         fprintf (fid, str1);
        %     else
        %         fid = fopen (pvd_filename, 'a');
        %     end
        % 
        %     if (fid < 0)
        %         error ('exportSolutionParaView: could not open file %s', pvd_filename);
        %     end
        % 
        %     ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
        %     if (isempty (ind)); ind = 1; end
        % 
        %     for iPatch = 1:obj.NumberPatches
        %         filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(currentIndex-1), '_', num2str (iPatch));
        %         filename_patch = cat (2, filename, '_', num2str(currentIndex-1), '_', num2str (iPatch));
        %         vts_pts = [];
        % 
        %         vts_pts(1, :, :) = obj.Geometry(iPatch).PlotPointsQuiver{1,1};
        %         vts_pts(2, :, :) = obj.Geometry(iPatch).PlotPointsQuiver{2,1};
        %         Az_patch = obj.MagneticPotential(obj.Spaces.gnum{iPatch}) / obj.Length;
        % 
        %         Bx = sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(2, :, :, :).* reshape(Az_patch, 1, 1, 1, []), 4);
        %         By = -sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(1, :, :, :).* reshape(Az_patch, 1, 1, 1, []), 4);
        %         angle = 0;
        %         Bxtemp = Bx*cos(angle) - By*sin(angle);
        %         Bytemp = By*cos(angle) + Bx*sin(angle);
        %         B = [Bxtemp; Bytemp];
        % 
        %         msh_to_vtk_vec (vts_pts, B, filename_patch, 'B');
        %         fprintf (fid, str2, dataList(currentIndex), iPatch, filename_patch_without_path);
        %     end
        % 
        %     if currentIndex == numel(dataList)
        %         fprintf (fid, str3);
        %     end
        %     fclose (fid);
        % end

        % set references to objects that hold the data?
    end
end
