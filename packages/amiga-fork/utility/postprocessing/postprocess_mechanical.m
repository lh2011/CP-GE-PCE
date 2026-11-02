function values = postprocess_mechanical(evalSpace, u, type)
    arguments
        evalSpace  
        u (:,1) double
        type {mustBeMember(type, ["Displacement", "DisplacementX", "DisplacementY", "Strain"])}
    end

    ndof_comp = prod (evalSpace.ndof_dir, 2);
    uDim{1} = u(1:ndof_comp(1));
    uDim{2} = u(ndof_comp(1)+(1:ndof_comp(2)));

    for iDim = 1:2        
        uEl{iDim} = zeros(size(evalSpace.scalar_spaces{iDim}.connectivity));
        nel{iDim} = size(evalSpace.scalar_spaces{iDim}.connectivity, 2);
        for iel = 1:nel{iDim}
            uEl{iDim}(:, iel) = uDim{iDim}(evalSpace.scalar_spaces{iDim}.connectivity(:, iel));
        end
    end

    switch type
        case "Displacement"
            ux = reshape(sum(evalSpace.scalar_spaces{1}.shape_functions.*permute(uEl{1}, [3,1,2]), 2), [], 1);
            uy = reshape(sum(evalSpace.scalar_spaces{2}.shape_functions.*permute(uEl{2}, [3,1,2]), 2), [], 1);
            values = sqrt(ux.^2 + uy.^2);
        case "DisplacementX"
            values = reshape(sum(evalSpace.scalar_spaces{1}.shape_functions.*permute(uEl{1}, [3,1,2]), 2), [], 1);
        case "DisplacementY"
            values = reshape(sum(evalSpace.scalar_spaces{2}.shape_functions.*permute(uEl{2}, [3,1,2]), 2), [], 1);
        case "Strain"
            nel = size(evalSpace.scalar_spaces{1}.shape_functions, 3);
            values = zeros(3, nel);
            % eps_xx
            values(1,:) = sum(evalSpace.scalar_spaces{1}.shape_function_gradients(1, :, :, :).*permute(uEl{1}, [3,4,1,2]), 3);
            % eps_yy
            values(2,:) = sum(evalSpace.scalar_spaces{2}.shape_function_gradients(2, :, :, :).*permute(uEl{2}, [3,4,1,2]), 3);
            % 2*epx_xy
            values(3,:) = sum(evalSpace.scalar_spaces{1}.shape_function_gradients(2, :, :, :).*permute(uEl{1}, [3,4,1,2]), 3) + ...
                          sum(evalSpace.scalar_spaces{2}.shape_function_gradients(1, :, :, :).*permute(uEl{2}, [3,4,1,2]), 3);
    end
end