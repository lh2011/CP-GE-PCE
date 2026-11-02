function varargout = int_D_DP_gradu_nu_gradv_v(spu, spv, spg, msh, u, v, dCdP, mat, len)
    nParams = size(dCdP, 3);
    uel = zeros(size(spu.connectivity));
    vel = zeros(size(spv.connectivity));
    dCdPel = zeros([size(spg.connectivity, 1), msh.ndim, size(dCdP, 3), msh.nel]);
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
        vel(:, iel) = v(spu.connectivity(:, iel));
        dCdPel(:, :, :, iel) = dCdP(spg.connectivity(:, iel), 1:msh.ndim, :);
    end

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4])/len;

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    gradu_u = sum(gradu.*permute(uel, [3,4,5,1,6,2]), 4);
    gradv_v = sum(gradv.*permute(vel, [3,4,1,5,6,2]), 3);
    B_mag = (sum(gradu_u.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;
    nu = mat.getReluctivity(B_mag);
    dnudB = mat.getNuPrime(B_mag);

    dK1dC = -weights.*jacdet.*nu.*sum(gradG.*gradv_v, 1).*gradu;
    dK2dC = -weights.*jacdet.*nu.*sum(gradG.*gradu, 1).*gradv_v;

    dK3dC = weights.*jacdet.*nu.*gradG.*sum(gradu.*gradv_v, 1);

    dK4dC = -weights.*jacdet.*sum(gradu.*gradv_v, 1).*dnudB./B_mag...
        .*gradu_u.*sum((gradG.*gradu_u), 1);

    dKdC = sum(dK1dC + dK2dC + dK3dC + dK4dC, 2);

    dKdP = squeeze(sum(dKdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5]));

    rows = repmat(permute(spu.connectivity, [3,1,2]), nParams, 1, 1);
    cols = repmat((1:nParams)', 1, spv.nsh_max, msh.nel);

    dKdP = reshape(dKdP, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, dKdP, spv.ndof, nParams);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dKdP;
    else
        error ('int_D_DP_gradu_nu_gradv_v: wrong number of output arguments')
    end
end