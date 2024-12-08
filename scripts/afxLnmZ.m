function afxLnmZ(image,fnameOut)
    if exist(fnameOut,'file')
        fprintf('afxLnmZ: %s already exists.\n', fnameOut);
        return;
    end

    [pth,~,~] = fileparts(fnameOut);
    if ~exist(pth,'dir'), mkdir(pth); end
    
    Vi = spm_vol(char(image));
    [X,XYZmm] = spm_read_vols(Vi);
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    clear Vtmp;

    gmMask = afxLoadMasks([],fullfile('masks','gmmask_20.nii'),XYZmm);
    
    Vo = struct(...
        'fname',   fnameOut,...
        'dim',     Vi.dim(1:3),...
        'dt',      [spm_type('int16') spm_platform('bigend')],...
        'pinfo',   [Inf Inf Inf]',...
        'mat',     Vi.mat,...
        'n',       1,...
        'descrip', 'afxMassUnivariate');
    
    Y = (X-mean(X(gmMask)))./std(X(gmMask));
    spm_write_vol(Vo,Y);
end