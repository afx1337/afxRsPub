function afxLnmTTest2(images1,images2,u,fnameOut)
%     if exist(fnameOut,'file')
%         fprintf('afxLnmZ: %s already exists.\n', fnameOut);
%         return;
%     end
    
    Vi1 = spm_vol(char(images1));
    Vi2 = spm_vol(char(images2));
  
    Vo = struct(...
        'fname',   fnameOut,...
        'dim',     Vi1(1).dim(1:3),...
        'dt',      [spm_type('int16') spm_platform('bigend')],...
        'pinfo',   [Inf Inf Inf]',...
        'mat',     Vi1(1).mat,...
        'n',       1,...
        'descrip', 'afxMassUnivariate');
    
    n1 = length(Vi1); % number of images
    n2 = length(Vi2); % number of images
    
    parfor i = 1:n1
        tmp = spm_read_vols(Vi1(i));
        X1(i,:) = tmp(:);
    end
    parfor i = 1:n2
        tmp = spm_read_vols(Vi2(i));
        X2(i,:) = tmp(:);
    end
    
    %implicit zero masking
    X1(X1==0) = NaN;
    X2(X2==0) = NaN;
    
    mask = ~(sum(isnan(X1),1) == size(X1,1) | sum(isnan(X2),1) == size(X2,1));
    fprintf('Two-sample t-test: analyzing %i voxels in %i vs %i images ...',length(mask),n1,n2);
    Y = nan(Vo.dim(1:3));
    [h,~,~,stats] = ttest2(X1(:,mask),X2(:,mask),u,'both');
    Y(mask) = stats.tstat.*h;
    spm_write_vol(Vo,Y);
    fprintf(' done\n');

end