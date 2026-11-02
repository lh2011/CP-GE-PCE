
function dbdP = int_D_DP_z_gradv_Br_bot(spv, spgP, msh, z, dCdP, material)
    nParams = size(dCdP, 3);
    zel = zeros(size(spv.connectivity));
    dCdPel = zeros([size(spgP.connectivity, 1), msh.ndim, nParams, msh.nel]);
    for iel = 1:msh.nel
        zel(:, iel) = z(spv.connectivity(:, iel));
        dCdPel(:, :, :, iel) = dCdP(spgP.connectivity(:, iel), 1:msh.ndim, :);
    end
    nu = material.getReluctivity();
    alpha = material.Angle(msh.geo_map(1,:,:), msh.geo_map(2,:,:));
    nbot = [-sin(alpha); cos(alpha)];
    nbot = permute(nbot, [1,2,4,5,6,3]);
    Br_bot = material.Br * nbot;
    
    gradG = permute(spgP.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradv_times_z = sum(gradv.*permute(zel, [3,4,1,5,6,2]), 3);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    db1dC = -weights.*jacdet.*nu.*sum(gradG.*Br_bot, 1).*gradv_times_z;
    db2dC = weights.*jacdet.*nu.*gradG.*sum(Br_bot.*gradv_times_z, 1);

    dbdC = sum(db1dC + db2dC, 2);

    dbdP = squeeze(sum(dbdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5,6]));
end
