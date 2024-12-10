function afxLnmTTest(images,fnameOut,nCond)
    if exist(fnameOut,'file')
        fprintf('afxLnmTTest: %s already exists.\n', fnameOut);
        return;
    end
    
    % read images
    Vi = spm_vol(char(images));
    % prepare output image
    Vo = struct(...
        'fname',   [],...
        'dim',     Vi(1).dim(1:3),...
        'dt',      [spm_type('int16') spm_platform('bigend')],...
        'pinfo',   [Inf Inf Inf]',...
        'mat',     Vi(1).mat,...
        'n',       1,...
        'descrip', 'afxMassUnivariate');
    n = length(Vi); % number of images
    
    % mean across conditions
    for i = 1:nCond:n
        for j = 1:nCond
            tmp = spm_read_vols(Vi(i+j-1));
            tmp2(j,:) = tmp(:);
        end
        X((i+nCond-1)/nCond,:) = mean(tmp2);
    end
    
    % implicit zero masking
    X(X==0) = NaN;
    mask = ~(sum(isnan(X),1) == size(X,1));
    
    fprintf('One-sample t-test: analyzing %i voxels in %i images ...',length(mask),size(X,1));
    
    % make output directories
    [pth,~,~] = fileparts(fnameOut);
    if ~exist(pth,'dir'), mkdir(pth); end
    
    % perform t-test (Y > 0?)
    [~,~,~,stats] = ttest(X(:,mask),0,[],'right');
    Y = nan(Vo.dim(1:3));
    Y(mask) = stats.tstat;
    Vo.fname = fnameOut;
    spm_write_vol(Vo,Y);

    fprintf(' done\n');
end