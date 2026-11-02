function varargout = int_D_DY_DC_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, v, mat, len)

    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end
    uel = permute(uel, [3,4,5,6,7,1,8,9,2]);
    vel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        vel(:, iel) = v(spu.connectivity(:, iel));
    end
    vel = permute(vel, [3,4,5,6,1,7,8,9,2]);
    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4]);
    gradNivi = sum(gradNi.*vel, 5);
    gradNj = permute(spu.shape_function_gradients, [1,5,6,2,7,3,8,9,4])/len;
    gradGk = permute(spg.shape_function_gradients, [1,5,6,2,7,8,3,9,4]);

    DkdgradNi = gradGk.*permute(gradNi, [2,1,3,4,5,6,7,8,9]);
    DkdgradNivi = sum(DkdgradNi.*vel, 5);
    DkdgradNj = gradGk.*permute(gradNj, [2,1,3,4,5,6,7,8,9]);
    DkdgradNjuj = sum(DkdgradNj.*uel, 6);

    TrDkd = permute(gradGk, [2,1,3,4,5,6,7,8,9]);


    gradNjuj = sum(gradNj.*uel, 6);
    B_mag = (sum(gradNjuj.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;
    nu = mat.getReluctivity(B_mag);
    dnudB = mat.getNuPrime(B_mag);
    dnu2dB2 = mat.getd2nudB2(B_mag);

    dBdCkd = - sum(DkdgradNjuj.*gradNjuj, 1)./B_mag;

    dBduj = sum(gradNjuj.*gradNj, 1)./B_mag;


    d2BdujdCkd = -sum(DkdgradNjuj.*gradNj, 1)./B_mag + ...
                -sum(gradNjuj.*DkdgradNj, 1)./B_mag +...
                - dBdCkd.*sum(gradNjuj.*gradNj, 1)./(B_mag).^2;

    J11 = sum(-weights.*dnudB.*dBduj.*sum(DkdgradNivi.*gradNjuj, 1).*jacdet, 4);
    J12 = sum(-weights.*dnudB.*dBduj.*sum(gradNivi.*DkdgradNjuj, 1).*jacdet, 4);
    J13 = sum(weights.*dnudB.*dBduj.*sum(gradNivi.*gradNjuj, 1).*TrDkd.*jacdet, 4);
    J14 = sum(weights.*dnu2dB2.*dBdCkd.*dBduj.*sum(gradNivi.*gradNjuj, 1).*jacdet, 4);
    J15 = sum(weights.*dnudB.*d2BdujdCkd.*sum(gradNivi.*gradNjuj, 1).*jacdet, 4);

    values = J11 + J12 + J13 + J14 + J15;

    nd = msh.ndim;
    ne =  1; %msh.ndim;
    ni = 1; %size(spv.connectivity, 1);
    nj = size(spu.connectivity, 1);
    nk = size(spg.connectivity, 1);
    nl = 1; %size(spg.connectivity, 1);
    nel = msh.nel;

    values = permute(values, [5,6,7,2,8,3,9,1,4]);

    % dimi = repmat(permute(spv.connectivity, [1,3,4,5,6,7,2]), 1 , nj, nk, nd, nl, ne, 1);
    dimj = repmat(permute(spu.connectivity, [3,1,4,5,6,7,2]), ni, 1 , nk, nd, nl, ne, 1);
    dimk = repmat(permute(spg.connectivity, [3,4,1,5,6,7,2]), ni, nj, 1 , nd, nl, ne, 1);
    dimd = repmat(permute((1:msh.ndim)',    [2,3,4,1]),       ni, nj, nk, 1 , nl, ne, nel);
    % diml = repmat(permute(spg.connectivity, [3,4,5,6,1,7,2]), ni, nj, nk, nd, 1 , ne, 1);
    % dime = repmat(permute((1:msh.ndim)',    [2,3,4,5,6,1]),   ni, nj, nk, nd, nl, 1 , nel);

    values = reshape(values, [], 1);
    dimd = reshape(dimd, [], 1);
    % dime = reshape(dime, [], 1);
    % dimi = reshape(dimi, [], 1);
    dimj = reshape(dimj, [], 1);
    dimk = reshape(dimk, [], 1);
    % diml = reshape(diml, [], 1);


    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([dimj, dimk, dimd], values, [spu.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [dimj, dimk, dimd];
        varargout{2} = values;
    else
        error ('int_D_DY_DC_gradu_u_nu_gradv_v: wrong number of output arguments')
    end
end
