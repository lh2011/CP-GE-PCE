function varargout = int_D_DC_DC_gradu_gradv(spu, spv, spg, msh)
    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4]);
    gradNj = permute(spu.shape_function_gradients, [1,5,6,2,7,3,8,9,4]);
    gradGk = permute(spg.shape_function_gradients, [1,5,6,2,7,8,3,9,4]);
    gradGl = permute(spg.shape_function_gradients, [1,5,6,2,7,8,9,3,4]);

    DkdDlegradNi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNi, [3,2,1,4,5,6,7,8,9]);
    DleDkdgradNi = permute(DkdDlegradNi, [1,3,2,4,5,6,8,7,9]);

    DkdDlegradNj = permute(DkdDlegradNi, [1,2,3,4,6,5,7,8,9]);
    DleDkdgradNj = permute(DleDkdgradNi, [1,2,3,4,6,5,7,8,9]);

    DkdgradNi = gradGk.*permute(gradNi, [2,1,3,4,5,6,7,8,9]);
    DlegradNi = permute(DkdgradNi, [1,3,2,4,5,6,8,7,9]);
    DkdgradNj = permute(DkdgradNi, [1,2,3,4,6,5,7,8,9]);
    DlegradNj = permute(DlegradNi, [1,2,3,4,6,5,7,8,9]);

    TrDkd = permute(gradGk, [2,1,3,4,5,6,7,8,9]);
    TrDle = permute(gradGl, [3,2,1,4,5,6,7,8,9]);

    TrDkdDle = TrDkd.*TrDle;
    TrDkdDle21 = TrDkdDle(:,2,1,:,:,:,:,:,:); % Temp value
    TrDkdDle(:,2,1,:,:,:,:,:,:) = TrDkdDle(:,1,2,:,:,:,:,:,:);  % Transpose d and e dimensions
    TrDkdDle(:,1,2,:,:,:,:,:,:) = TrDkdDle21;


    K11 = sum(weights.*sum(DleDkdgradNi.*gradNj, 1).*jacdet, 4);
    K12 = sum(weights.*sum(DkdDlegradNi.*gradNj, 1).*jacdet, 4);
    K13 = sum(weights.*sum(DkdgradNi.*DlegradNj, 1).*jacdet, 4);
    K14 = sum(-weights.*sum(DkdgradNi.*gradNj, 1).*jacdet.*TrDle, 4);

    dKdCdC1 = K11 + K12 + K13 + K14;

    dKdCdC2 = permute(dKdCdC1, [1,2,3,4,6,5,7,8,9]);

    K31 = sum(-weights.*sum(DlegradNi.*gradNj, 1).*jacdet.*TrDkd, 4);
    K32 = sum(-weights.*sum(gradNi.*DlegradNj, 1).*jacdet.*TrDkd, 4);
    K33 = sum(weights.*sum(gradNi.*gradNj, 1).*jacdet.*TrDle.*TrDkd, 4);
    K34 = sum(-weights.*sum(gradNi.*gradNj, 1).*jacdet.*TrDkdDle, 4);

    dKdCdC3 = K31 + K32 + K33 + K34;

    
    values = dKdCdC1 + dKdCdC2 + dKdCdC3;

    nd = msh.ndim;
    ne = msh.ndim;
    ni = size(spv.connectivity, 1);
    nj = size(spu.connectivity, 1);
    nk = size(spg.connectivity, 1);
    nl = size(spg.connectivity, 1);
    nel = msh.nel;

    values = permute(values, [5,6,7,2,8,3,9,1,4]);

    dimi = repmat(permute(spv.connectivity, [1,3,4,5,6,7,2]), 1 , nj, nk, nd, nl, ne, 1);
    dimj = repmat(permute(spu.connectivity, [3,1,4,5,6,7,2]), ni, 1 , nk, nd, nl, ne, 1);
    dimk = repmat(permute(spg.connectivity, [3,4,1,5,6,7,2]), ni, nj, 1 , nd, nl, ne, 1);
    dimd = repmat(permute((1:msh.ndim)',    [2,3,4,1]),       ni, nj, nk, 1 , nl, ne, nel);
    diml = repmat(permute(spg.connectivity, [3,4,5,6,1,7,2]), ni, nj, nk, nd, 1 , ne, 1);
    dime = repmat(permute((1:msh.ndim)',    [2,3,4,5,6,1]),   ni, nj, nk, nd, nl, 1 , nel);

    values = reshape(values, [], 1);
    dimd = reshape(dimd, [], 1);
    dime = reshape(dime, [], 1);
    dimi = reshape(dimi, [], 1);
    dimj = reshape(dimj, [], 1);
    dimk = reshape(dimk, [], 1);
    diml = reshape(diml, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([dimi, dimj, dimk, dimd, diml, dime], values, [spv.ndof, spu.ndof, spg.ndof, msh.ndim, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [dimi, dimj, dimk, dimd, diml, dime];
        varargout{2} = values;
    else
        error ('int_D_DC_DC_gradu_gradv: wrong number of output arguments')
    end
end
