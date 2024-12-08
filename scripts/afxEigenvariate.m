function Y = afxEigenvariate(y)
    % taken from spm_regions.m
    [m,n]   = size(y);
    if m > n
        [~,s,v] = svd(y'*y);
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
    else
        [~,s,u] = svd(y*y');
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
    end
    d       = sign(sum(v));
    u       = u*d;
    %v       = v*d;
    Y       = u*sqrt(s(1)/n);
end