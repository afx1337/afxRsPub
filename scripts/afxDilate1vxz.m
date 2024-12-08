function afxDilate1vxz()
    [t,~] = spm_select([1 Inf],'image','Select masks ...');
    for i = 1:size(t,1)
        fname = t(i,:);
        V = spm_vol(fname);
        Y = spm_read_vols(V);
        Y(:,:,1:end-1) = Y(:,:,1:end-1) | Y(:,:,2:end);
        Y(:,:,2:end) = Y(:,:,2:end) | Y(:,:,1:end-1);
        [p,f,~] = fileparts(fname);
        pNew = fullfile(p,'dilatedZ1vx');
        if ~exist(pNew,'dir')
            mkdir(pNew);
        end
        fnameNew = fullfile(pNew,strcat('dilZ1vx_',f,'.nii'));
        V.fname = fnameNew;
        spm_write_vol(V,Y);
    end
end