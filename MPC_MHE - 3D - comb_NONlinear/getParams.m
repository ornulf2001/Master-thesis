function data = getParams(parameters)
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end