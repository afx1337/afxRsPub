function Y = afxEigenvariate(y)
    % Y = afxEigenvariate(y)
    % First eigenvariate, sign-aligned to ROI mean
    % y: [T x V]

    % remove voxel-wise mean for PCA
    y = y - mean(y,1);

    % ROI mean (reference)
    m = mean(y,2);

    if size(y,2) > 230 && size(y,2) < 2000
        Y = afxEigenvariate_new(y);
    else
        Y = afxEigenvariate_old(y);
    end

    % enforce mean-consistent sign
    if dot(Y, m) < 0
        Y = -Y;
    end
end

function Y = afxEigenvariate_new(y)
    % first singular component
    try
        [u,s,~] = svds(y,1);
        s = s(1,1);
    catch
        [u,S,~] = svd(y,'econ');
        u = u(:,1);
        s = S(1,1);
    end

    % scale like SPM
    Y = u * (s / sqrt(size(y,2)));
end

function Y = afxEigenvariate_old(y)
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
        %v       = y'*u/sqrt(s(1));
    end
    %d       = sign(sum(v));
    %u       = u*d;
    %%v       = v*d;
    Y       = u*sqrt(s(1)/n);
end
