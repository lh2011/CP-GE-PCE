function varargout = int_D_DC_DC_gradv_v_Br_bot(spv, spg, msh, material, z)
    zel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        zel(:, iel) = z(spv.connectivity(:, iel));
    end
    zel = permute(zel, [3,4,5,6,1,7,8,9,2]);

    nu = material.getReluctivity();
    alpha = material.Angle(0,0);
    Br_bot = nu*material.Br * [-sin(alpha); cos(alpha)];

    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4]);
    gradNizi = sum(gradNi.*zel, 5);
    gradGk = permute(spg.shape_function_gradients, [1,5,6,2,7,8,3,9,4]);
    gradGl = permute(spg.shape_function_gradients, [1,5,6,2,7,8,9,3,4]);

    % DkdDlegradNi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNi, [3,2,1,4,5,6,7,8,9]);
    % DleDkdgradNi = permute(DkdDlegradNi, [1,3,2,4,5,6,8,7,9]);
    DkdDlegradNizi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNizi, [3,2,1,4,5,6,7,8,9]);
    DleDkdgradNizi = gradGl.*permute(gradGk, [3,2,1,4,5,6,7,8,9]).*permute(gradNizi, [2,1,3,4,5,6,7,8,9]);

    % DkdgradNi = gradGk.*permute(gradNi, [2,1,3,4,5,6,7,8,9]);
    % DlegradNi = permute(DkdgradNi, [1,3,2,4,5,6,8,7,9]);
    DkdgradNizi = gradGk.*permute(gradNizi, [2,1,3,4,5,6,7,8,9]);
    DlegradNizi = gradGl.*permute(gradNizi, [3,2,1,4,5,6,7,8,9]);

    TrDkd = permute(gradGk, [2,1,3,4,5,6,7,8,9]);
    TrDle = permute(gradGl, [3,2,1,4,5,6,7,8,9]);

    TrDkdDle = TrDkd.*TrDle;
    TrDkdDle21 = TrDkdDle(:,2,1,:,:,:,:,:,:); % Temp value
    TrDkdDle(:,2,1,:,:,:,:,:,:) = TrDkdDle(:,1,2,:,:,:,:,:,:);  % Transpose d and e dimensions
    TrDkdDle(:,1,2,:,:,:,:,:,:) = TrDkdDle21;


    b11 = sum(weights.*sum(DleDkdgradNizi.*Br_bot, 1).*jacdet, 4);
    b12 = sum(weights.*sum(DkdDlegradNizi.*Br_bot, 1).*jacdet, 4);
    b13 = sum(-weights.*sum(DkdgradNizi.*Br_bot, 1).*jacdet.*TrDle, 4);

    dbdCdC1 = b11 + b12 + b13;

    b21 = sum(-weights.*sum(DlegradNizi.*Br_bot, 1).*jacdet.*TrDkd, 4);
    b22 = sum(weights.*sum(gradNizi.*Br_bot, 1).*jacdet.*TrDle.*TrDkd, 4);
    b23 = sum(-weights.*sum(gradNizi.*Br_bot, 1).*jacdet.*TrDkdDle, 4);

    dbdCdC2 = b21 + b22 + b23;

    values = dbdCdC1 + dbdCdC2;

    nd = msh.ndim;
    ne = msh.ndim;
    ni = 1; %size(spv.connectivity, 1);
    nj = 1; %size(spu.connectivity, 1);
    nk = size(spg.connectivity, 1);
    nl = size(spg.connectivity, 1);
    nel = msh.nel;

    values = permute(values, [5,6,7,2,8,3,9,1,4]);

    % dimi = repmat(permute(spv.connectivity, [1,3,4,5,6,7,2]), 1 , nj, nk, nd, nl, ne, 1);
    % dimj = repmat(permute(spu.connectivity, [3,1,4,5,6,7,2]), ni, 1 , nk, nd, nl, ne, 1);
    dimk = repmat(permute(spg.connectivity, [3,4,1,5,6,7,2]), ni, nj, 1 , nd, nl, ne, 1);
    dimd = repmat(permute((1:msh.ndim)',    [2,3,4,1]),       ni, nj, nk, 1 , nl, ne, nel);
    diml = repmat(permute(spg.connectivity, [3,4,5,6,1,7,2]), ni, nj, nk, nd, 1 , ne, 1);
    dime = repmat(permute((1:msh.ndim)',    [2,3,4,5,6,1]),   ni, nj, nk, nd, nl, 1 , nel);

    values = reshape(values, [], 1);
    dimd = reshape(dimd, [], 1);
    dime = reshape(dime, [], 1);
    % dimi = reshape(dimi, [], 1);
    % dimj = reshape(dimj, [], 1);
    dimk = reshape(dimk, [], 1);
    diml = reshape(diml, [], 1);


    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([dimk, dimd, diml, dime], values, [spg.ndof, msh.ndim, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [dimk, dimd, diml, dime];
        varargout{2} = values;
    else
        error ('int_D_DC_DC_gradv_Br_bot_times_z: wrong number of output arguments')
    end
end
