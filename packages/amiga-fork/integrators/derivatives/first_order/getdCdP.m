function dCdP = getdCdP(geometryFunction, allParameters, testParameters, spaceGeo, Steps)
    paramNames = fieldnames(testParameters);
    nParameters = numel(fieldnames(testParameters));
    if ~exist("Steps", "var")
        Steps = ones(nParameters, 1)*1e-6;
    end

    C0 = getC(geometryFunction, allParameters, testParameters, spaceGeo);
    dCdP = zeros([size(C0), nParameters]);
    
    for iParam = 1:nParameters
        Step = Steps(iParam);
        testParametersPerturbed = testParameters;
        testParametersPerturbed.(paramNames{iParam}) = testParametersPerturbed.(paramNames{iParam})+Step;
        C1 = getC(geometryFunction, allParameters, testParametersPerturbed, spaceGeo);

        dCdP(:, :, iParam) = (C1 - C0)/(Step);

        % testParametersPerturbed2 = testParameters;
        % testParametersPerturbed2.(paramNames{iParam}) = testParametersPerturbed2.(paramNames{iParam})-Step;
        % C2 = getC(geometryFunction, allParameters, testParametersPerturbed2, spaceGeo);
        % 
        % dCdP(:, :, iParam) = (C1 - C2)/(2*Step);
    end
end
