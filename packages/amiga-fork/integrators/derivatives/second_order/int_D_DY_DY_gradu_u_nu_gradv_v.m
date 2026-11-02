function varargout = int_D_DY_DY_gradu_u_nu_gradv_v(spu, spv, msh, u, v, mat, len)
    uel = zeros(size(spu.connectivity));
    vel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
        vel(:, iel) = v(spu.connectivity(:, iel));
    end
    uel = permute(uel, [3,4,5,6,7,1,8,9,2]);
    vel = permute(vel, [3,4,5,6,1,7,8,9,2]);
    %  1  2 3 4 5 6 7 8 9
    % sum d e q i j k l el
    % Definitions:
    jacdet = permute(msh.jacdet, [3,4,5,1,6,7,8,9,2]);
    weights = permute(msh.quad_weights, [3,4,5,1,6,7,8,9,2]);

    gradNi = permute(spv.shape_function_gradients, [1,5,6,2,3,7,8,9,4])/len;
    gradNj = permute(spu.shape_function_gradients, [1,5,6,2,7,3,8,9,4])/len;

    gradNjuj = sum(gradNj.*uel, 6);
    gradNivi = sum(gradNi.*vel, 5);
    B_mag = (sum(gradNjuj.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;

    dnudB = mat.getNuPrime(B_mag);
    d2nudB2 = mat.getd2nudB2(B_mag);

    dBduj = sum(gradNjuj.*gradNj, 1)./B_mag;
    dBdui = sum(gradNjuj.*gradNi, 1)./B_mag;

    dBduiduj = sum(gradNi.*gradNj, 1)./B_mag - sum(gradNjuj.*gradNi, 1)./(B_mag.^2).* dBduj;

    J11 = sum(len.*weights.*d2nudB2.*dBduj.*dBdui.*sum(gradNjuj.*gradNivi, 1).*jacdet, 4);
    J12 = sum(len.*weights.*dnudB.*dBduiduj.*sum(gradNjuj.*gradNivi, 1).*jacdet, 4);

    J13 = sum(len.*weights.*dnudB.*dBdui.*sum(gradNivi.*gradNj, 1).*jacdet, 4) + ...
          sum(len.*weights.*dnudB.*dBduj.*sum(gradNivi.*gradNi, 1).*jacdet, 4);

    values = J11 + J12 + J13;

    nd = 1;
    ne =  1; %msh.ndim;
    ni = size(spu.connectivity, 1);
    nj = size(spu.connectivity, 1);
    nk = 1;
    nl = 1; %size(spg.connectivity, 1);
    nel = msh.nel;

    values = permute(values, [5,6,7,2,8,3,9,1,4]);

    dimi = repmat(permute(spu.connectivity, [1,3,4,5,6,7,2]), 1 , nj, nk, nd, nl, ne, 1);
    dimj = repmat(permute(spu.connectivity, [3,1,4,5,6,7,2]), ni, 1 , nk, nd, nl, ne, 1);

    values = reshape(values, [], 1);

    dimi = reshape(dimi, [], 1);
    dimj = reshape(dimj, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (dimi, dimj, values, spu.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = dimi;
        varargout{2} = dimj;
        varargout{3} = values;
    else
        error ('int_D_DY_DY_gradu_u_nu_gradv_v: wrong number of output arguments')
    end
end
