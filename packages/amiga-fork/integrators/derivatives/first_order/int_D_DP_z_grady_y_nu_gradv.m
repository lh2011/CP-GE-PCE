
function zdKdP = int_D_DP_z_grady_y_nu_gradv(spu, spv, spg, msh, y, z, dCdP, mat)
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

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    gradv_times_z = sum(gradv.*permute(zel, [3,4,1,5,6,2]), 3);

    gradu_times_u = sum(gradu.*permute(yel, [3,4,5,1,6,2]), 4);
    B_mag = (sum(gradu_times_u.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;

    nu = mat.getReluctivity(B_mag);
    dnudB = mat.getNuPrime(B_mag);

    zdK1dC = -weights.*jacdet.*nu.*sum(gradG.*gradv_times_z, 1).*gradu_times_u;
    zdK2dC = -weights.*jacdet.*nu.*sum(gradG.*gradu_times_u, 1).*gradv_times_z;

    zdK3dC = weights.*jacdet.*nu.*gradG.*sum(gradu_times_u.*gradv_times_z, 1);

    zdK4dC = -weights.*jacdet.*sum(gradu_times_u.*gradv_times_z, 1).*dnudB./B_mag...
        .*gradu_times_u.*sum((gradG.*gradu_times_u), 1); 

    zdKdC = sum(zdK1dC + zdK2dC + zdK3dC + zdK4dC, 2);

    zdKdP = squeeze(sum(zdKdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5,6]));
end