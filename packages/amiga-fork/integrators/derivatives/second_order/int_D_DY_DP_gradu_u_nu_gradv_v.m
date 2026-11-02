function varargout = int_D_DY_DP_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, v, dCdP, mat, len)

    dKdYdC = int_D_DY_DC_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, v, mat, len);
    
    dKdYdP = double(ttt(dKdYdC, sptensor(dCdP(:,1:msh.ndim,:)), [2,3], [1,2]));
    
   
    if (nargout == 1 || nargout == 0)
        varargout{1} = dKdYdP;
    elseif (nargout == 3)
        [rows, cols, vals] = find(dKdYdP);
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = vals;
    else
        error ('int_D_DY_DP_gradu_u_nu_gradv_v: wrong number of output arguments')
    end
end
