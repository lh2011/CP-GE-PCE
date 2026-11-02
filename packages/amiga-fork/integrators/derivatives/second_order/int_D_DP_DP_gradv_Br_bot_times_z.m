function varargout = int_D_DP_DP_gradv_Br_bot_times_z(spv, spg, msh, z, dCdP, d2CdP2, mat)

    % First order derivative
    dbdC = int_D_DC_gradv_v_Br_bot(spv, spg, msh, mat, z);
    % Second order derivative
    dbdCdC = int_D_DC_DC_gradv_v_Br_bot(spv, spg, msh, mat, z);
    
    % Product rule
    % dbdPdP1 = ttt(dbdC, sptensor(d2CdP2(:,1:msh.ndim,:,:)), [2,3], [1,2]);
    dbdPdP1 = squeeze(sum(full(dbdC).*d2CdP2(:,1:msh.ndim,:,:), [1,2]));

    dbdPdP2 = ttt(sptensor(ttt(dbdCdC, sptensor(dCdP(:,1:msh.ndim,:)), [3,4], [1,2])), ...
        sptensor(dCdP(:,1:msh.ndim,:)), [1,2], [1,2]);

    dbdPdP = dbdPdP1 + double(dbdPdP2);

    if (nargout == 1 || nargout == 0)
        varargout{1} = dbdPdP;
    elseif (nargout == 3)
        [rows, cols, values] = find(dbdPdP);
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error ('int_D_DP_DP_gradv_Br_bot_times_z: wrong number of output arguments')
    end
end
