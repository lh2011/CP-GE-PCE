function d2CdP2 = getd2CdP2(geometryFunction, allParameters, testParameters, spaceGeo, Steps)
    paramNames = fieldnames(testParameters);
    nParameters = numel(fieldnames(testParameters));
    if ~exist("Steps", "var")
        Steps = ones(nParameters, 1)*1e-6;
    end

    dCdP0 = getdCdP(geometryFunction, allParameters, testParameters, spaceGeo, Steps);
    d2CdP2 = zeros([size(dCdP0), nParameters]);
    

    for iParam = 1:nParameters
        Step = Steps(iParam);
        testParametersPerturbed = testParameters;
        testParametersPerturbed.(paramNames{iParam}) = testParametersPerturbed.(paramNames{iParam})+Step;
        dCdP1 = getdCdP(geometryFunction, allParameters, testParametersPerturbed, spaceGeo, Steps);

        d2CdP2(:,:,:,iParam) = (dCdP1-dCdP0)/(Step);

        % testParametersPerturbed = testParameters;
        % testParametersPerturbed.(paramNames{iParam}) = testParametersPerturbed.(paramNames{iParam})-Step;
        % dCdP2 = getdCdP(geometryFunction, allParameters, testParametersPerturbed, spaceGeo, Steps);
        % d2CdP2(:,:,:,iParam) = (dCdP1-dCdP2)/(2*Step);
    end
end