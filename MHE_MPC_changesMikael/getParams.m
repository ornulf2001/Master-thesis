function data = getParams()
    persistent params

    if isempty(params)
        parameters;
    end
    data = params;
end