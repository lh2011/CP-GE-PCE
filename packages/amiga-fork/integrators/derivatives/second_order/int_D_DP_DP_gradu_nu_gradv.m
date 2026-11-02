function varargout = int_D_DP_DP_gradu_nu_gradv(spu, spv, spg, msh, dCdP, d2CdP2, u, mat)

    % First order derivative
    dKdC = int_D_DC_gradu_nu_gradv(spu, spv, spg, msh, u, mat);
    % Second order derivative
    dKdCdC = int_D_DC_DC_gradu_nu_gradv(spu, spv, spg, msh, u, mat);
    
    % Product rule
    dKdPdP1 = ttt(dKdC, sptensor(d2CdP2(:,1:msh.ndim,:,:)), [3,4], [1,2]);

    dKdPdP2 = ttt(sptensor(ttt(dKdCdC, sptensor(dCdP(:,1:msh.ndim,:)), [5,6], [1,2])), ...
        sptensor(dCdP(:,1:msh.ndim,:)), [3,4], [1,2]);

    dKdPdP = dKdPdP1 + dKdPdP2;

    if (nargout == 1 || nargout == 0)
        varargout{1} = dKdPdP;
    elseif (nargout == 2)
        [indices, values] = find(dKdPdP);
        varargout{1} = indices;
        varargout{2} = values;
    else
        error ('int_D_DP_DP_gradu_nu_gradv: wrong number of output arguments')
    end
end
