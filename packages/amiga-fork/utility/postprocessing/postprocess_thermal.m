function values = postprocess_thermal(evalSpace, a, type)
    arguments
        evalSpace  
        a (:,1) double
        type {mustBeMember(type, ["T"])}
    end
    aEl = zeros(size(evalSpace.connectivity));
    nel = size(evalSpace.connectivity, 2);
    for iel = 1:nel
        aEl(:, iel) = a(evalSpace.connectivity(:, iel));
    end

    switch type
        case "T"
            values = reshape(sum(evalSpace.shape_functions.*permute(aEl, [3,1,2]), 2), [], 1);
    end

end