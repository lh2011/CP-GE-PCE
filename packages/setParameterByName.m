function parameters = setParameterByName(parameters, parameterNames, name, value)
    % find index for parameter
    idx = find(parameterNames == string(name));
    assert(~isempty(idx), 'Optimization parameter "%s" was not found.', name);
    assert(isscalar(idx), 'Optimization parameter "%s" occurs more than once.', name);
    % assign value
    parameters(idx) = value;
end
