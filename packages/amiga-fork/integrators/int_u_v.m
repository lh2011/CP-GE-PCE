% tbd, own implementation

function varargout = int_u_v(spu, spv, msh, func)
    if ~exist("func", "var")
        func = @(x,y,z) ones(size(x));
    end
    if isnumeric(func)
        func = @(x,y,z) func*ones(size(x));
    end
    for idim = 1:msh.rdim
      x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
    end
    coef = permute(func(x{:}), [3,1,4,5,2]);

    shpu = reshape (spu.shape_functions, spu.ncomp, msh.nqn, 1, spu.nsh_max, msh.nel); % [dim, nquad, v, u, nel]
    shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, 1, msh.nel); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    values = reshape(sum(sum(shpu.*shpv, 1).*jacdet.*coef, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('int_u_v: wrong number of output arguments')
    end
end
