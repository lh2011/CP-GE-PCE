
function zdKdP = int_D_DP_z_grady_y_gradv(spu, spv, spg, msh, y, z, dCdP)
    nParams = size(dCdP, 3);
    yel = zeros(size(spu.connectivity));
    zel = zeros(size(spu.connectivity));
    dCdPel = zeros([size(spg.connectivity, 1), msh.ndim, nParams, msh.nel]);
    for iel = 1:msh.nel
        yel(:, iel) = y(spu.connectivity(:, iel));
        zel(:, iel) = z(spu.connectivity(:, iel));
        dCdPel(:, :, :, iel) = dCdP(spg.connectivity(:, iel), 1:msh.ndim, :);
    end

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);

    gradv_times_z = sum(gradv.*permute(zel, [3,4,1,5,6,2]), 3);
    grady_times_y = sum(gradu.* permute(yel, [3,4,5,1,6,2]), 4);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    dK1dC = -weights.*jacdet.*sum(gradG.*gradv_times_z, 1).*grady_times_y;
    dK2dC = -weights.*jacdet.*sum(gradG.*grady_times_y, 1).*gradv_times_z;
    dK3dC = weights.*jacdet.*gradG.*sum(grady_times_y.*gradv_times_z, 1);

    zdKdC = sum(dK1dC + dK2dC + dK3dC, 2);

    zdKdP = squeeze(sum(zdKdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5,6]));
end
