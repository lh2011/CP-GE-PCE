function dVdP = int_D_DP_Omega(spg, msh, dCdP)
    
    dVdC = int_D_DC_Omega(spg, msh);

    dVdP = reshape(sum(full(dVdC).*dCdP(:,1:msh.ndim, :), [1,2]), [], 1);
end

