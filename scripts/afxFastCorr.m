function rho = afxFastCorr(x, y)
    % rho = afxFastCorr(x, y)
    %
    % x and y must be matrixes or vectors with same
    % number of elements. All elements must be real.
    % afxFastCorr performs no further checks.

    % vectorize
    x = x(:);
    y = y(:);

    % mean-center (in place)
    x = x - mean(x);
    y = y - mean(y);

    % pearson correlation coefficient
    rho = dot(x, y) / (vecnorm(x) .* vecnorm(y));
end
