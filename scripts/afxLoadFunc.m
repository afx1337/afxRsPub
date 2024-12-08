function [y,XYZmm,dim,mat] = afxLoadFunc(func)

    fprintf('   Loading functional data ... ')
    % read functional data and convert to one dimensional vector

    [~,~,e] = fileparts(func(1,:));
    if strcmp(e,'.mat')
        raw = load(func(1,:));
        dim = raw.dim;
        mat = raw.mat;
        [R,C,P] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
        RCP     = [R(:)';C(:)';P(:)';ones(1,numel(R))];
        XYZmm   = mat(1:3,:)*RCP;
        XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
        y = nan(size(raw.Y,1),length(raw.brainMask));
        y(:,raw.brainMask) = double(raw.Y);
    else
        % new (~4 seconds for a NKI data set)
        for  i = 1:size(func,1)
            tmp = nifti(func(i,:));
            if size(tmp.dat,4) > 1
                parfor j = 1:size(tmp.dat,4)
                    y(j,:) = reshape(tmp.dat(:,:,:,j),1,[]);
                end
                break
            else
                y(i,:) = tmp.dat(:);
            end
        end

        % load first functional image for MNI coordinates in mm and for
        % dimensions etc.
        Vfunc = spm_vol(func(1,:));
        [~,XYZmm] = spm_read_vols(spm_vol(Vfunc));
        XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
        dim = Vfunc.dim;
        mat = Vfunc.mat;
    end
    fprintf('done\n');
end