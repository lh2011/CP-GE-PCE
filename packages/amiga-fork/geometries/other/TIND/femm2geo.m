function femm2geo(femmName, geoName, args)
    arguments
        femmName 
        geoName 
        args.Scaling = 1;
        args.MeshSize = 1e-3;
        args.RefineArea = [0, 0];
        args.RefineRatio = 0.1;
    end
    femmdata = read_femmdata(num2str(femmName));
    points = femmdata.points*args.Scaling;
    segments = femmdata.segments;
    arcsegments = femmdata.arcsegments;

    fid = fopen(geoName, 'w');

    % write all points
    for iPoint = 1:size(points, 1)
        p = points(iPoint, 1:2);
        r = (p(1)^2 + p(2)^2)^0.5;
        if r > args.RefineArea(1) && r < args.RefineArea(2)
            ms = args.MeshSize * args.RefineRatio;
        else
            ms = args.MeshSize;
        end
        fprintf(fid, 'Point(%d) = {%.6f, %.6f, 0, %.6f};\n', iPoint, p(1), p(2), ms);
    end
    
    % write line segments
    for iSeg = 1:size(segments, 1)
        p1 = segments(iSeg, 1)+1;
        p2 = segments(iSeg, 2)+1;
        fprintf(fid, 'Line(%d) = {%d, %d};\n', iSeg, p1, p2);
    end

    % write circular segments
    nSeg = iSeg;
    for iArcSeg = 1:size(arcsegments, 1)
        p1 = points(arcsegments(iArcSeg,1)+1, 1:2);
        p2 = points(arcsegments(iArcSeg,2)+1, 1:2);
        alpha = arcsegments(iArcSeg, 3);
        
        mid = 0.5*(p1+p2);
        secant = p2 - p1;
        nArc = [-secant(2), secant(1)];
        nArc = nArc/vecnorm(nArc);
        dMid = vecnorm(secant);
        R = dMid/(2*sind(alpha/2));
        l = sqrt(R.^2 - (dMid/2).^2);
        center = mid + l*nArc; % Centre of circle point

        % write center points as new points
        iPoint = iPoint + 1;
        fprintf(fid, 'Point(%d) = {%.6f, %.6f, 0, %.6f};\n', iPoint, center(1), center(2), args.MeshSize);
        
        % split circle if there is a 180 degree segment
        if alpha+1e-12 > 180
            pmid = center - nArc*R;
            iPoint = iPoint + 1;
            % write half point of circle
            fprintf(fid, 'Point(%d) = {%.6f, %.6f, 0, %.6f};\n', iPoint, pmid(1), pmid(2), args.MeshSize);
            % write two circles
            nSeg = nSeg + 1;
            fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', nSeg, arcsegments(iArcSeg,1)+1, iPoint-1, iPoint);
            nSeg = nSeg + 1;
            fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', nSeg, iPoint, iPoint-1, arcsegments(iArcSeg,2)+1);
        else
            nSeg = nSeg + 1;
            % write circle
            fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', nSeg, arcsegments(iArcSeg,1)+1, iPoint, arcsegments(iArcSeg,2)+1);
        end
    end
    fclose(fid);
end