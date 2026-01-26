function afxLnmSecondlevelCalc(files, fMean, fT, nCond)
    % afxLnmSecondlevelCalc(files, fMean, fT)
    
    % load data
    [y,~,dim,mat] = afxLoadFunc(char(files));

    % conditions mitteln
    y = squeeze(mean(reshape(y,nCond,[],size(y,2)),1));
    
    % mean, sd and n
    mu = mean(y, 1);
    mu2 = mean(y.^2,1);
    n = size(y, 1);
    %sd = std(y, 0, 1);
    %t = mu ./ (sd ./ sqrt(n));
    var = (mu2 - mu.^2) * (n/(n-1));  % unbiased variance
    t = mu ./ sqrt(var / n);
    
    % write mean
    afxVolumeWrite(fMean,mu,dim,'int16',mat,'Lesion network map: mean',true);

    % calculate and write t-statistic
    afxVolumeWrite(fT,t,dim,'int16',mat,'Lesion network map: t-test',true);
end