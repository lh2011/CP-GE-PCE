% adjusted copy from nrbplot
function preview_nurbs(srf, regions)
    figureObj = findobj('type', 'figure', 'name', 'preview');
    if ishandle(figureObj)
        figure(figureObj);
    else
        figure('Name', 'preview', 'NumberTitle', 'off');
    end
    regionNames = fieldnames(regions);
    % check region definitions
    for iRegion = 1:numel(regionNames)
        for iRegion2 = 1:numel(regionNames)
            if iRegion ~= iRegion2
                if any(intersect(regions.(regionNames{iRegion}), regions.(regionNames{iRegion2})))
                    error("Regions are not uniquely defined!");
                end
            end
        end
    end

    clf
    hold on
    for iSrf = 1:numel(srf)
        for iReg = 1:numel(regionNames)
            if ismember(iSrf, regions.(regionNames{iReg}))
                break;
            end
            if iReg == numel(regionNames)
                iReg = 0;
            end
        end
        if iReg == 0
            color = 'black';
            warning("Region not specified for surface %i ", iSrf);
        else
            color = num2str(TUDa_getColor(iReg));
        end
        nrbplotcol(srf(iSrf), [10, 10], 'color', color);
    end
    axis equal
    view(2)

    function nrbplotcol (nurbs, subd, varargin)
    
        nargs = nargin;
        if nargs < 2
            if ~exist('subd','var')
                subd = nurbs.number.* nurbs.order;
            else
                error ('Need a NURBS to plot and the number of subdivisions!');
            end
        elseif rem(nargs+2,2)
            error ('Param value pairs expected')
        end
        
        if (any (nurbs.order < 2))
            warning ('The plot with order smaller than 2 may not be correct')
        end
        
        % Default values
        light='off';
        cmap='summer';
        plot_col='green';
        % Recover Param/Value pairs from argument list
        for i=1:2:nargs-2
            Param = varargin{i};
            Value = varargin{i+1};
            if (~ischar (Param))
                error ('Parameter must be a string')
            elseif size(Param,1)~=1
                error ('Parameter must be a non-empty single row string.')
            end
            switch lower (Param)
                case 'light'
                    light = lower (Value);
                    if (~ischar (light))
                        error ('light must be a string.')
                    elseif ~(strcmp(light,'off') | strcmp(light,'on'))
                        error ('light must be off | on')
                    end
                case 'colormap'
                    if ischar (Value)
                        cmap = lower(Value);
                    elseif size (Value, 2) ~= 3
                        error ('colormap must be a string or have exactly three columns.')
                    else
                        cmap=Value;
                    end
                case 'color'
                    if ischar (Value)
                        plot_col = lower(Value);
                    elseif size (Value, 2) ~= 3
                        error ('colormap must be a string or have exactly three columns.')
                    else
                        plot_col=Value;
                    end
                otherwise
                    error ('Unknown parameter: %s', Param)
            end
        end
        
        %colormap (cmap);
        
        % convert the number of subdivisions in number of points
        subd = subd+1;
        
        % plot the curve or surface
        if (iscell (nurbs.knots))
            if (size (nurbs.knots,2) == 2) % plot a NURBS surface
                knt = nurbs.knots;
                order = nurbs.order;
                p = nrbeval (nurbs, {linspace(knt{1}(order(1)),knt{1}(end-order(1)+1),subd(1)) ...
                    linspace(knt{2}(order(2)),knt{2}(end-order(2)+1),subd(2))});
                if (strcmp (light,'on'))
                    % light surface
                    surfl (squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)));
                    shading interp;
                else
                    surf (squeeze (p(1,:,:)), squeeze (p(2,:,:)), squeeze (p(3,:,:)),'FaceAlpha',0.75,'FaceColor',plot_col,'EdgeColor','none');
                    %shading faceted;
                end
            elseif (size (nurbs.knots,2) == 3) % plot the boundaries of a NURBS volume
                bnd = nrbextract (nurbs);
                hold_flag = ishold;
                nrbplotcol (bnd(1), subd(2:3), varargin{:});
                hold on
                nrbplotcol (bnd(2), subd(2:3), varargin{:});
                nrbplotcol (bnd(3), subd([1 3]), varargin{:});
                nrbplotcol (bnd(4), subd([1 3]), varargin{:});
                nrbplotcol (bnd(5), subd(1:2), varargin{:});
                nrbplotcol (bnd(6), subd(1:2), varargin{:});
        
                if (~hold_flag)
                    hold off
                end
        
            else
                error ('nrbplot: some argument is not correct')
            end
        else
            % plot a NURBS curve
            order = nurbs.order;
            p = nrbeval (nurbs, linspace (nurbs.knots(order), nurbs.knots(end-order+1), subd));
        
            if (any (nurbs.coefs(3,:)))
                % 3D curve
                plot3 (p(1,:), p(2,:), p(3,:));
                grid on;
            else
                % 2D curve
                plot (p(1,:), p(2,:), "Color", [0.5,0.5,0.5]);
            end
        end
        % axis equal;
    
    end
end

