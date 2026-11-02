function varargout = int_D_DC_DC_gradu_u_gradv_v(spu, spv, spg, msh, u, z)
    % slower version:
    % D_DC_DC = int_D_DC_DC_gradu_gradv(spu, spv, spg, msh, u);
    % 
    % D_DC_DC_times_u = ttv(D_DC_DC, u, 2);
    % 
    % if (nargout == 1 || nargout == 0)
    %     varargout{1} = D_DC_DC_times_u;
    % elseif (nargout == 2)
    %     [indices1, values1] = find(D_DC_DC_times_u);
    %     varargout{1} = indices1;
    %     varargout{2} = values1;
    % else
    %     error ('int_D_DC_DC_gradu_gradv_times_u: wrong number of output arguments')
    % end


    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end
    uel = permute(uel, [3,4,5,6,7,1,8,9,2]);

    zel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        zel(:, iel) = z(spv.connectivity(:, iel));
    end
    zel = permute(zel, [3,4,5,6,1,7,8,9,2]);
    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4]);
    gradNj = permute(spu.shape_function_gradients, [1,5,6,2,7,3,8,9,4]);
    gradGk = permute(spg.shape_function_gradients, [1,5,6,2,7,8,3,9,4]);
    gradGl = permute(spg.shape_function_gradients, [1,5,6,2,7,8,9,3,4]);

    gradNjuj = sum(gradNj.*uel, 6);
    gradNizi = sum(gradNi.*zel, 5);

    % DkdDlegradNi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNi, [3,2,1,4,5,6,7,8,9]);
    % DleDkdgradNi = permute(DkdDlegradNi, [1,3,2,4,5,6,8,7,9]);
    
    % DkdDlegradNizi = sum(DkdDlegradNi.*zel, 5);
    % DleDkdgradNizi = sum(DleDkdgradNi.*zel, 5);
    DkdDlegradNizi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNizi, [3,2,1,4,5,6,7,8,9]);
    DleDkdgradNizi = gradGl.*permute(gradGk, [3,2,1,4,5,6,7,8,9]).*permute(gradNizi, [2,1,3,4,5,6,7,8,9]);

    % DkdDlegradNj = permute(DkdDlegradNi, [1,2,3,4,6,5,7,8,9]);
    % DleDkdgradNj = permute(DleDkdgradNi, [1,2,3,4,6,5,7,8,9]);

    DkdDlegradNjuj = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNjuj, [3,2,1,4,5,6,7,8,9]);
    DleDkdgradNjuj = gradGl.*permute(gradGk, [3,2,1,4,5,6,7,8,9]).*permute(gradNjuj, [2,1,3,4,5,6,7,8,9]);
    % DkdDlegradNjuj = sum(DkdDlegradNj.*uel, 6);
    % DleDkdgradNjuj = sum(DleDkdgradNj.*uel, 6);

    % DkdgradNi = gradGk.*permute(gradNi, [2,1,3,4,5,6,7,8,9]);
    % DkdgradNizi = sum(DkdgradNi.*zel, 5);
    DkdgradNizi = gradGk.*permute(gradNizi,  [2,1,3,4,5,6,7,8,9]);

    % DlegradNi = permute(DkdgradNi, [1,3,2,4,5,6,8,7,9]);
    % DlegradNizi = sum(DlegradNi.*zel, 5);
    DlegradNizi = gradGl.*permute(gradNizi,  [3,2,1,4,5,6,7,8,9]);

    % DkdgradNj = permute(DkdgradNi, [1,2,3,4,6,5,7,8,9]);
    % DkdgradNjuj = sum(DkdgradNj.*uel, 6);
    DkdgradNjuj =gradGk.*permute(gradNjuj,  [2,1,3,4,5,6,7,8,9]);

    % DlegradNj = permute(DlegradNi, [1,2,3,4,6,5,7,8,9]);
    % DlegradNjuj = sum(permute(DlegradNi, [1,2,3,4,6,5,7,8,9]).*uel, 6);
    DlegradNjuj = gradGl.*permute(gradNjuj,  [3,2,1,4,5,6,7,8,9]);

    TrDkd = permute(gradGk, [2,1,3,4,5,6,7,8,9]);
    TrDle = permute(gradGl, [3,2,1,4,5,6,7,8,9]);

    TrDkdDle = TrDkd.*TrDle;
    TrDkdDle21 = TrDkdDle(:,2,1,:,:,:,:,:,:); % Temp value
    TrDkdDle(:,2,1,:,:,:,:,:,:) = TrDkdDle(:,1,2,:,:,:,:,:,:);  % Transpose d and e dimensions
    TrDkdDle(:,1,2,:,:,:,:,:,:) = TrDkdDle21;


    K11 = sum(weights.*sum(DleDkdgradNizi.*gradNjuj, 1).*jacdet, 4);
    K12 = sum(weights.*sum(DkdDlegradNizi.*gradNjuj, 1).*jacdet, 4);
    K13 = sum(weights.*sum(DkdgradNizi.*DlegradNjuj, 1).*jacdet, 4);
    K14 = sum(-weights.*sum(DkdgradNizi.*gradNjuj, 1).*jacdet.*TrDle, 4);

    dKdCdC1 = K11 + K12 + K13 + K14;

    K21 = sum(weights.*sum(DlegradNizi.*DkdgradNjuj, 1).*jacdet, 4);
    K22 = sum(weights.*sum(gradNizi.*DleDkdgradNjuj, 1).*jacdet, 4);
    K23 = sum(weights.*sum(gradNizi.*DkdDlegradNjuj, 1).*jacdet, 4);
    K24 = sum(-weights.*sum(gradNizi.*DkdgradNjuj, 1).*jacdet.*TrDle, 4);

    dKdCdC2 = K21 + K22 + K23 + K24;

    K31 = sum(-weights.*sum(DlegradNizi.*gradNjuj, 1).*jacdet.*TrDkd, 4);
    K32 = sum(-weights.*sum(gradNizi.*DlegradNjuj, 1).*jacdet.*TrDkd, 4);
    K33 = sum(weights.*sum(gradNizi.*gradNjuj, 1).*jacdet.*TrDle.*TrDkd, 4);
    K34 = sum(-weights.*sum(gradNizi.*gradNjuj, 1).*jacdet.*TrDkdDle, 4);

    dKdCdC3 = K31 + K32 + K33 + K34;

    values = dKdCdC1 + dKdCdC2 + dKdCdC3;

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
        error ('int_D_DC_DC_gradu_u_gradv_v: wrong number of output arguments')
    end
end
