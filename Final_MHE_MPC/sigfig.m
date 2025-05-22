function x_rounded = sigfig(x, n)
    % Rounds every element in x to n significant figures
    x_rounded = zeros(size(x));
    for idx = 1:numel(x)
        if x(idx) == 0
            x_rounded(idx) = 0;
        else
            d = n - floor(log10(abs(x(idx)))) - 1;
            x_rounded(idx) = round(x(idx), d);
        end
    end
end