function values = postprocess_magnetic(evalSpace, a, type)
    arguments
        evalSpace  
        a (:,1) double
        type {mustBeMember(type, ["Az", "Bx", "By", "Bmag"])}
    end
    aEl = zeros(size(evalSpace.connectivity));
    nel = size(evalSpace.connectivity, 2);
    for iel = 1:nel
        aEl(:, iel) = a(evalSpace.connectivity(:, iel));
    end

    switch type
        case "Az"
            values = reshape(sum(evalSpace.shape_functions.*permute(aEl, [3,1,2]), 2), [], 1);
        case "Bx"
            values = reshape(sum(evalSpace.shape_function_gradients(2,:,:,:).*permute(aEl, [3,4,1,2]), 3), [], 1);
        case "By"
            values = -reshape(sum(evalSpace.shape_function_gradients(1,:,:,:).*permute(aEl, [3,4,1,2]), 3), [], 1);
        case "Bmag"
            values = reshape(sum(sum(evalSpace.shape_function_gradients.*permute(aEl, [3,4,1,2]), 3).^2, 1).^0.5, [], 1);
    end

end