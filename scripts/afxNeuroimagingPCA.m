function afxNeuroimagingPCA(files,k,outdir)
    % load data
	[y,~,dim,mat] = afxLoadFunc(char(files));
    
    % generate implizit mask and discard all data outside the mask
    mask = mean(abs(y),1) > eps;
    y = y(:,mask);
    
    % pca
    [U, explained, mu, V, explainedVoxelPC] = afxFastPCA(y, k);
    clear y;
    
    % output dir
    if ~exist(outdir,'dir'), mkdir(outdir); end
    
    % save loadings
    fac = zeros(1,k);
    for i = 1:k
        % scale loadings to std
        %fac(i) = std(U(:,i)); % Varianzverschiebung von U nach V
        fac(i) = 1/std(V(:,i)); % ~ konstant
        V(:,i) = V(:,i).*fac(i);
        afxVolumeWrite(fullfile(outdir,sprintf('V_%03d.nii',i)),afxUnmask(V(:,i),mask),dim,'int16',mat,'PCA Loadings',true);
        afxVolumeWrite(fullfile(outdir,sprintf('expl_%03d.nii',i)),afxUnmask(explainedVoxelPC(:,i),mask),dim,'int16',mat,'Voxelwise fraction of variance explained by component',true);
    end
    
    % scale U to match scaled loadings
    Ufac = U./fac;
    
    % save weights, factors and explained variance
    save(fullfile(outdir,'result.mat'),'U','explained','fac','Ufac');
    
    % save mean
    afxVolumeWrite(fullfile(outdir,'mu.nii'),afxUnmask(mu,mask),dim,'int16',mat,'Mean',true);
end

function y_full = afxUnmask(y,mask)
    y_full = zeros(1, numel(mask),'like',y); 
    y_full(:,mask) = y;
end