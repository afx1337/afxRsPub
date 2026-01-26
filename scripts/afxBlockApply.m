function X = afxBlockApply(X, fun, blockSize)
    % afxBlockApply  Wendet blockweise eine Funktion auf Spalten von X an
    %
    %   Y = afxBlockApply(X, fun, blockSize)
    %
    %   X         : T x V Matrix (double)
    %   fun       : Function handle, z. B. @(x) filtfilt(b,a,x)
    %   blockSize : Anzahl der Spalten pro Block (default 1000)
    %
    %   Y         : T x V Ergebnismatrix

    if nargin < 3 || isempty(blockSize)
        blockSize = 1000;
    end

    [~, cols] = size(X);

    for i = 1:blockSize:cols
        idx = i:min(i+blockSize, cols);
        X(:,idx) = fun(X(:,idx));
    end
end
