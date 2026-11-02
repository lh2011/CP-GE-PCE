function varargout = int_f_v (spv, msh, func)
    if ~exist("func", "var")
        func = @(x,y,z) repmat(reshape(ones(size(x)), [1, size(x)]), spv.ncomp, 1, 1);
    end
    if isnumeric(func)
        func = @(x,y,z) func*repmat(reshape(ones(size(x)), [1, size(x)]), spv.ncomp, 1, 1);
    end
    for idim = 1:msh.rdim
        x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
    end
    coeff = func(x{:});  % [ncomp, nqn, nel]
    coeff = reshape (coeff, spv.ncomp, msh.nqn, 1, msh.nel);
    
    shpv  = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel); % [ncomp, nqn, nv, nel]
    
    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,2]);
    
    values = reshape(sum(sum(shpv.*coeff, 1) .* jacdet, 2), [], 1);
    
    rows = reshape(spv.connectivity, [], 1);
    cols = ones(size(rows));
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, 1);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('int_f_v: wrong number of output arguments')
    end
end

