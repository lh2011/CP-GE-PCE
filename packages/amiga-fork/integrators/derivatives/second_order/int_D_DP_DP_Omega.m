
function dVdPdP = int_D_DP_DP_Omega(spg, msh, dCdP, d2CdP2)

    nParams = size(dCdP, 3);
    dCdPel = zeros([size(spg.connectivity, 1), msh.ndim, nParams, msh.nel]);
    for iel = 1:msh.nel
        dCdPel(:, :, :, iel) = dCdP(spg.connectivity(:, iel), 1:msh.ndim, :);
    end

    gradGk = permute(spg.shape_function_gradients, [4,2,3,1]);
    gradGl = permute(spg.shape_function_gradients, [4,2,5,6,3,1]);

    jacdet = permute(msh.jacdet, [2,1]);
    weights = permute(msh.quad_weights, [2,1]);

    d2VdC2_1 = sum(weights.*jacdet.*gradGl.*gradGk, 2);

    traceGlGk = gradGk.*gradGl;
    traceGlGk21 = traceGlGk(:,:,:,2,:,1); % temp value
    traceGlGk(:,:,:,2,:,1) = traceGlGk(:,:,:,1,:,2);
    traceGlGk(:,:,:,1,:,2) = traceGlGk21;

    d2VdC2_2 = -sum(weights.*jacdet.*traceGlGk, 2);

    % Second order derivative
    d2VdC2 = squeeze(d2VdC2_1 + d2VdC2_2);

    % First order derivative
    dVdC = int_D_DC_Omega(spg, msh);
    
    %% Product rule
    dVdPdP1 = squeeze(sum(sum(d2VdC2.*permute(dCdPel, [4,1,2,5,6,3]), [2,3]) .*permute(dCdPel, [4,5,6,1,2,7,3]), [1,4,5]));
    dVdPdP2 = squeeze(sum(full(dVdC).*d2CdP2(:,1:2,:,:), [1,2]));
    
    dVdPdP = dVdPdP1 + dVdPdP2;
end
