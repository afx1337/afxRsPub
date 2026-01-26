function [y,excl,rmsFd] = afxDenoising(y,subjectMasks,brainMask,rpFile,options)

    % if >options.unsmoothed< is true, tissue signal will be extracted from
    % unsmoothed data (y2), confound removal and filtering will be performed
    % also for unsmoothed data

    fprintf('   Denoising ...')

    % calculate rp and rp'
    rp = dlmread(rpFile);
    rpp = [zeros(1,size(rp,2)); diff(rp)];
    motion = nan(size(rp,1),0);
    if options.regressRp
        motion = [ rp rp.^2 rpp rpp.^2 ];
    end
    
    % calculate tissue signal regressors
    ts = nan(size(y,1),0);
    nanMask = any(isnan(y),1);
    for i = 1:length(options.regressTs)
        % generate mask
        combinedProb = subjectMasks(:,1:3)*options.regressTs(i).tpm';
        tissueMask = combinedProb > options.regressTs(i).thresh;
        tissueMask(nanMask) = false; % to be safe from voxels outside the brain in the HCP data set
        if nnz(tissueMask) == 0
            error('Empty tissue mask for tissue signal regression/residualization.')
        end
        
        % get tissue nuisance regressors
        if options.regressTs(i).pca > 0
            % compute principle component analysis
            tmp = afxFastPCA(y(:,tissueMask), options.regressTs(i).pca);
        else
            % compute weighted mean (in a memory efficient way)
            tmp = weightedMeanBlockwise(y, combinedProb, tissueMask, 1000);
        end
        ts = [ts tmp];
    end
    clear combinedProb tissueMask nanMask;
   
    % generate counfound matrix, mean center and norm confounds
    confounds = [ motion ts ];
    confounds = confounds - mean(confounds,1);
    confounds = confounds ./ std(confounds,0,1);
    % confound removal via GLM/multiple regression: y = X*beta + epsilon
    % optimized for avoiding RAM peaks
    if size(confounds,2) > 1
        X = single([confounds ones(size(confounds,1),1)]);
        [Q,~] = qr(X,0);     
        y(:,brainMask) = afxBlockApply(y(:,brainMask), @(Yblock) Yblock - Q * (Q' * Yblock), 1000);
    end

    % apply filter
    if options.filter(1) ~= 0 || ~isinf(options.filter(2))
        y(:,brainMask) = afxFilter(options.TR,options.filter,y(:,brainMask));
    end
    
    % motion scrubbing
    % calculate framewise displacment (FD)
    % (maximum displacment of a voxel within a 50 mm sphere)
    FD = sum(abs([rpp(:,1:3) rpp(:,4:6)*50]),2);
    % exclude scans with more than 0.5 mm FD
    % see Power et al. (2012) and Varkruti et al. (2016)
    extMov = FD > options.threshFD;
    y = y(~extMov,:);

    if nargout > 2, excl = nnz(extMov); end
    if nargout > 3, rmsFd = sqrt(mean(FD.^2)); end
    
    fprintf(' done\n')
    fprintf('      > %i of %i scans were discarded due to motion scrubbing\n',nnz(extMov),numel(FD));
end


function tmp = weightedMeanBlockwise(y, weights, mask, blockSize)
    % Maskierte Indizes
    idxMask = find(mask);
    nMask   = numel(idxMask);

    % Initialisierung
    N   = size(y,1);
    tmp = zeros(N,1,'like',y);
    wsum = zeros(1,1,'like',weights);

    % Blockweise Akkumulation
    for k = 1:blockSize:nMask
        idx = idxMask(k:min(k+blockSize-1, nMask));
        w   = weights(idx);

        tmp  = tmp + y(:,idx) * w(:);
        wsum = wsum + sum(w,'native');
    end

    % Normierung
    tmp = tmp / wsum;
end
