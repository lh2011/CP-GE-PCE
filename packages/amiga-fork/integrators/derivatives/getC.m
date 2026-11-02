function ControlPoints = getC(geometryFunction, allParameters, testParameters, spaceGeo)
    paramNames = fieldnames(testParameters);
    for iParam = 1:numel(paramNames)
        allParameters.(paramNames{iParam}) = testParameters.(paramNames{iParam});
    end
    unpackedParams = namedargs2cell(allParameters);
    srf = geometryFunction(unpackedParams{:});
    
    ControlPoints = zeros(spaceGeo.ndof, 4);
    for iPatch  = 1:spaceGeo.npatch
        ind_loc = spaceGeo.gnum{iPatch};
        ControlPoints(ind_loc, 1) = reshape(srf(iPatch).coefs(1, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 2) = reshape(srf(iPatch).coefs(2, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 3) = reshape(srf(iPatch).coefs(3, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 4) = reshape(srf(iPatch).coefs(4, :, :, :), [], 1);
    end
end
