% adjusted copy of sp_eval to return precomputed spaces and meshes

function [mshEval, spaceEval] = sp_eval_pnts (space, geometry, npts)

    ndim = numel (space.knots);
    
    endpoints = zeros (2, ndim);
    if (isfield (geometry, 'nurbs'))
        nurbs = geometry.nurbs;
        if (ndim == 1)
            nurbs.knots = {nurbs.knots};
        end
        for idim=1:ndim
            endpoints(:,idim) = nurbs.knots{idim}([nurbs.order(idim), end-nurbs.order(idim)+1]);
        end
        clear nurbs
    elseif (isfield (struct(space), 'knots'))
        degree = space.degree;
        for idim=1:ndim
            endpoints(:,idim) = space.knots{idim}([degree(idim)+1, end-degree(idim)]);
        end
    else
        endpoints(2,:) = 1;
    end
    
    if (iscell (npts))
        pts = npts;
        npts = cellfun (@numel, pts);
        degree = space.degree;
        for idim = 1:ndim
            knt_aux = space.knots{idim}(degree(idim)+1:end-degree(idim));
            ind = find (pts{idim}(1)>=knt_aux(1:end-1), 1, 'last');
            endpoints(1,idim) = knt_aux(ind);
            ind = find (pts{idim}(end)>=knt_aux(1:end-1), 1, 'last');
            endpoints(2,idim) = knt_aux(ind+1);
        end
    elseif (isvector (npts))
        if (numel (npts) == 1)
            npts = npts * ones (1,ndim);
        end
        for idim = 1:ndim
            pts{idim} = linspace (endpoints(1,idim), endpoints(2,idim), npts(idim));
        end
    end
    
    for jj = 1:ndim
        pts{jj} = pts{jj}(:)';
        if (numel (pts{jj}) > 1)
            brk{jj} = [endpoints(1,jj), pts{jj}(1:end-1) + diff(pts{jj})/2, endpoints(2,jj)];
        else
            brk{jj} = endpoints(:,jj).';
        end
    end
    
    msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
    mshEval = msh_precompute(msh);
    sp  = space.constructor (msh);
    spaceEval = sp_precompute(sp, mshEval, 'value', true, 'gradient', true);
end
