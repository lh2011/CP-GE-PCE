
function varargout = int_D_DC_DC_Omega(spg, msh)

    gradGk = permute(spg.shape_function_gradients, [4,2,3,1]);
    gradGl = permute(spg.shape_function_gradients, [4,2,5,6,3,1]);

    jacdet = permute(msh.jacdet, [2,1]);
    weights = permute(msh.quad_weights, [2,1]);

    dV1dC = sum(weights.*jacdet.*gradGl.*gradGk, 2);

    traceGlGk = gradGk.*gradGl;
    traceGlGk21 = traceGlGk(:,:,:,2,:,1); % temp value
    traceGlGk(:,:,:,2,:,1) = traceGlGk(:,:,:,1,:,2);
    traceGlGk(:,:,:,1,:,2) = traceGlGk21;

    dV2dC = -sum(weights.*jacdet.*traceGlGk, 2);

    dVdC = squeeze(dV1dC + dV2dC);

    rows = repmat(permute(spg.connectivity, [2,1]), 1, 1, msh.ndim, size(spg.connectivity, 1), msh.ndim);
    dims1 = repmat(permute((1:msh.ndim)', [2,3,1]), msh.nel, size(spg.connectivity, 1), 1, size(spg.connectivity, 1), msh.ndim);

    cols = repmat(permute(spg.connectivity, [2,3,4,1]), 1, size(spg.connectivity, 1), msh.ndim, 1, msh.ndim);
    dims2 = repmat(permute((1:msh.ndim)', [2,3,4,5,1]), msh.nel, size(spg.connectivity, 1), msh.ndim, size(spg.connectivity, 1), 1);

    dVdC = reshape(dVdC, [], 1);
    rows = reshape(rows, [], 1);
    dims1 = reshape(dims1, [], 1);
    cols = reshape(cols, [], 1);
    dims2 = reshape(dims2, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows, dims1, cols, dims2], dVdC, [spg.ndof, msh.ndim, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows, dims1, cols, dims2];
        varargout{2} = dVdC;
    else
        error ('int_D_DC_DC_Omega: wrong number of output arguments')
    end
end
