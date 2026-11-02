function varargout = int_D_DY_DC_gradu_u_nu_gradv(spu, spv, spg, msh, u, mat)

    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end
    uel = permute(uel, [3,4,5,6,7,1,8,9,2]);
    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4]);
    gradNj = permute(spu.shape_function_gradients, [1,5,6,2,7,3,8,9,4]);
    gradGk = permute(spg.shape_function_gradients, [1,5,6,2,7,8,3,9,4]);
    % gradGl = permute(spg.shape_function_gradients, [1,5,6,2,7,8,9,3,4]);

    % DkdDlegradNi = gradGk.*permute(gradGl, [2,1,3,4,5,6,7,8,9]).*permute(gradNi, [3,2,1,4,5,6,7,8,9]);
    % DleDkdgradNi = permute(DkdDlegradNi, [1,3,2,4,5,6,8,7,9]);
    % 
    % DkdDlegradNj = permute(DkdDlegradNi, [1,2,3,4,6,5,7,8,9]);
    % DkdDlegradNjuj = sum(DkdDlegradNj.*uel, 6);
    % DleDkdgradNj = permute(DleDkdgradNi, [1,2,3,4,6,5,7,8,9]);
    % DleDkdgradNjuj = sum(DleDkdgradNj.*uel, 6);

    DkdgradNi = gradGk.*permute(gradNi, [2,1,3,4,5,6,7,8,9]);
    DlegradNi = permute(DkdgradNi, [1,3,2,4,5,6,8,7,9]);
    DkdgradNj = permute(DkdgradNi, [1,2,3,4,6,5,7,8,9]);
    DkdgradNjuj = sum(DkdgradNj.*uel, 6);
    DlegradNj = permute(DlegradNi, [1,2,3,4,6,5,7,8,9]);
    DlegradNjuj = sum(permute(DlegradNi, [1,2,3,4,6,5,7,8,9]).*uel, 6);

    TrDkd = permute(gradGk, [2,1,3,4,5,6,7,8,9]);
    % TrDle = permute(gradGl, [3,2,1,4,5,6,7,8,9]);

    % TrDkdDle = TrDkd.*TrDle;
    % TrDkdDle21 = TrDkdDle(:,2,1,:,:,:,:,:,:); % Temp value
    % TrDkdDle(:,2,1,:,:,:,:,:,:) = TrDkdDle(:,1,2,:,:,:,:,:,:);  % Transpose d and e dimensions
    % TrDkdDle(:,1,2,:,:,:,:,:,:) = TrDkdDle21;

    gradNjuj = sum(gradNj.*uel, 6);
    B_mag = (sum(gradNjuj.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;
    nu = mat.getNuNonlinear(B_mag);
    dnudB = mat.getNuPrimeNonlinear(B_mag);
    dnu2dB2 = mat.getd2nudB2(B_mag);

    dBdCkd = -sum(gradNjuj.*gradGk, 1).*permute(gradNjuj, [2,1,3,4,5,6,7,8,9])./B_mag;

    dBduj = sum(gradNjuj.*gradNj, 1)./B_mag;

    % dBdCle = permute(dBdCkd, [1,3,2,4,5,6,8,7,9]);
    % 
    % dBdCkddCle = 1./B_mag.*(sum(sum(DlegradNj.*uel, 6).*sum(DkdgradNj.*uel, 6), 1) ...
    %                        + sum(gradNjuj.*sum(DleDkdgradNj.*uel, 6), 1) ...
    %                        + sum(gradNjuj.*sum(DkdDlegradNj.*uel, 6), 1) ...
    %                        - dBdCkd.*dBdCle);

    d2BdudCkd = -sum(DkdgradNjuj.*gradNj, 1)./B_mag + ...
                -sum(gradNjuj.*DkdgradNj, 1)./B_mag +...
                - dBdCkd.*sum(gradNjuj.*gradNj, 1)./(B_mag).^2;

    J11 = sum(-weights.*dnudB.*dBduj.*sum(DkdgradNi.*gradNjuj, 1).*jacdet, 4);
    J12 = sum(-weights.*dnudB.*dBduj.*sum(gradNi.*DkdgradNjuj, 1).*jacdet, 4);
    J13 = sum(weights.*dnudB.*dBduj.*sum(gradNi.*gradNjuj, 1).*TrDkd.*jacdet, 4);
    J14 = sum(weights.*dnu2dB2.*dBdCkd.*dBduj.*sum(gradNi.*gradNjuj, 1).*jacdet, 4);
    J15 = sum(weights.*dnudB.*d2BdudCkd.*sum(gradNi.*gradNjuj, 1).*jacdet, 4);

    values = J11 + J12 + J13 + J14 + J15;

    nd = msh.ndim;
    ne =  1; %msh.ndim;
    ni = size(spv.connectivity, 1);
    nj = size(spu.connectivity, 1);
    nk = size(spg.connectivity, 1);
    nl = 1; %size(spg.connectivity, 1);
    nel = msh.nel;

    values = permute(values, [5,6,7,2,8,3,9,1,4]);

    dimi = repmat(permute(spv.connectivity, [1,3,4,5,6,7,2]), 1 , nj, nk, nd, nl, ne, 1);
    dimj = repmat(permute(spu.connectivity, [3,1,4,5,6,7,2]), ni, 1 , nk, nd, nl, ne, 1);
    dimk = repmat(permute(spg.connectivity, [3,4,1,5,6,7,2]), ni, nj, 1 , nd, nl, ne, 1);
    dimd = repmat(permute((1:msh.ndim)',    [2,3,4,1]),       ni, nj, nk, 1 , nl, ne, nel);
    % diml = repmat(permute(spg.connectivity, [3,4,5,6,1,7,2]), ni, nj, nk, nd, 1 , ne, 1);
    % dime = repmat(permute((1:msh.ndim)',    [2,3,4,5,6,1]),   ni, nj, nk, nd, nl, 1 , nel);

    values = reshape(values, [], 1);
    dimd = reshape(dimd, [], 1);
    % dime = reshape(dime, [], 1);
    dimi = reshape(dimi, [], 1);
    dimj = reshape(dimj, [], 1);
    dimk = reshape(dimk, [], 1);
    % diml = reshape(diml, [], 1);


    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([dimi, dimj, dimk, dimd], values, [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [dimi, dimj, dimk, dimd];
        varargout{2} = values;
    else
        error ('int_D_DY_DC_gradu_u_nu_gradv: wrong number of output arguments')
    end
end
