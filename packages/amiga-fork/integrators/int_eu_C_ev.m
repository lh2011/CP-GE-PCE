% tbd, own implementation, only for 2D

function varargout = int_eu_C_ev(spu, spv, msh, C)
    
    gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], 1, msh.nqn, 1, spu.nsh_max, msh.nel);
    gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], 1, msh.nqn, spv.nsh_max, 1, msh.nel);

    jacdet = permute(msh.jacdet .* msh.quad_weights, [1,3,4,2]);

    Bu = zeros(3, 2, 1, msh.nqn, 1, spu.nsh_max, msh.nel);
    Bu(1,1,:,:,:,:,:) = gradu(1,1,:,:,:,:,:);
    Bu(2,2,:,:,:,:,:) = gradu(2,2,:,:,:,:,:);
    Bu(3,1,:,:,:,:,:) = gradu(1,2,:,:,:,:,:);
    Bu(3,2,:,:,:,:,:) = gradu(2,1,:,:,:,:,:);

    Bv = zeros(3, 2, 1, msh.nqn, spu.nsh_max, 1, msh.nel);
    Bv(1,1,:,:,:,:,:) = gradv(1,1,:,:,:,:,:);
    Bv(2,2,:,:,:,:,:) = gradv(2,2,:,:,:,:,:);
    Bv(3,1,:,:,:,:,:) = gradv(1,2,:,:,:,:,:);
    Bv(3,2,:,:,:,:,:) = gradv(2,1,:,:,:,:,:);
    
    % Bv'*C*Bu
    BvCBu = sum(permute(Bv, [2,3,1,4,5,6,7]).*sum(Bu.*permute(C, [2,3,1]), 1), 3);
    % many values in BvCBu are zero, due to stupid ordering of multi-dimensinoal spaces in GeoPDEs.
    BvCBu = reshape(sum(BvCBu, [1,2]), msh.nqn, spv.nsh_max, spu.nsh_max, msh.nel);

    values = reshape(sum(BvCBu.*jacdet, 1), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('int_eu_C_ev: wrong number of output arguments')
    end
end
