function varargout = int_D_DY_gradu_u_nu_gradv(spu, spv, msh, u, material)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    gradu = permute(spu.shape_function_gradients, [1,2,5,3,4]); % [dim, nquad, v, u, nel]
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    gradNi_ui = sum(gradu.*permute(uel, [3,4,5,1,2]), 4);

    B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;

    dnu_dB = material.getNuPrimeNonlinear(B_mag);

    values = reshape(sum(sum(gradv.*gradNi_ui, 1).*(sum(gradNi_ui.*gradu, 1)).*jacdet.*dnu_dB./B_mag, 2), [], 1);
    % values = reshape(sum(sum(gradv.*gradu, 1).*nu.*jacdet + sum(gradv.*gradNi_ui, 1).*(sum(gradNi_ui.*gradu, 1)).*jacdet.*dnu_dB./B_mag, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error ('int_D_DY_gradu_u_nu_gradv: wrong number of output arguments')
    end
end


