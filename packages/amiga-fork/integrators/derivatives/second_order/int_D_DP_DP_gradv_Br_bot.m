function varargout = int_D_DP_DP_gradv_Br_bot(spv, spg, msh, dCdP, d2CdP2, mat)

    % First order derivative
    dbdC = int_D_DC_gradv_Br_bot(spv, spg, msh, mat);
    % Second order derivative
    dbdCdC = int_D_DC_DC_gradv_Br_bot(spv, spg, msh, mat);
    
    % Product rule
    dbdPdP1 = ttt(dbdC, sptensor(d2CdP2(:,1:msh.ndim,:,:)), [2,3], [1,2]);

    dbdPdP2 = ttt(sptensor(ttt(dbdCdC, sptensor(dCdP(:,1:msh.ndim,:)), [4,5], [1,2])), ...
        sptensor(dCdP(:,1:msh.ndim,:)), [2,3], [1,2]);

    dbdPdP = dbdPdP1 + dbdPdP2;

    if (nargout == 1 || nargout == 0)
        varargout{1} = dbdPdP;
    elseif (nargout == 2)
        [indices, values] = find(dbdPdP);
        varargout{1} = indices;
        varargout{2} = values;
    else
        error ('int_D_DP_DP_gradv_Br_bot: wrong number of output arguments')
    end
end
