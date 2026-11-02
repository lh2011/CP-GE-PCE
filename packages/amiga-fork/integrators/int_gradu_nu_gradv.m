% TBD 

function varargout = int_gradu_nu_gradv(spu, spv, msh, a, material)
    ael = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        ael(:, iel) = a(spu.connectivity(:, iel));
    end

    gradu = permute(spu.shape_function_gradients, [1,2,5,3,4]); % [dim, nquad, v, u, nel]
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    gradNi_ui = sum(gradu.*permute(ael, [3,4,5,1,2]), 4);

    B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
    nu = material.getReluctivity(B_mag);

    values = reshape(sum(sum(gradu.*gradv, 1).*jacdet.*nu, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('int_gradu_nu_gradv: wrong number of output arguments')
    end
end
