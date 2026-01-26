function [y,XYZmm,dim,mat] = afxLoadFunc(func, varargin)
    % Load NIfTI images with optional caching and dynamic cutoff
    %
    % Optional parameters:
    %   'precision' : 'single' (default) or 'logical'

    p = inputParser;
    addParameter(p,'precision','single',@(x) ismember(x,{'single','logical'}));
    parse(p,varargin{:});
    precision = p.Results.precision;

    fprintf('   Loading functional data ... ')

    [~,~,e] = fileparts(func(1,:));
    if strcmp(e,'.mat')
        raw = load(func(1,:));
        dim = raw.dim;
        mat = raw.mat;
        [R,C,P] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
        RCP     = [R(:)';C(:)';P(:)';ones(1,numel(R))];
        XYZmm   = mat(1:3,:)*RCP;
        XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
        y = nan(size(raw.Y,1),length(raw.brainMask),precision);
        if strcmp(precision,'logical')
            y(:,raw.brainMask) = raw.Y~=0;
        else
            y(:,raw.brainMask) = single(raw.Y);
        end
    else
        %% --- Estimate memory requirement ---
        tmp = nifti(func(1,:));
        nVolumes = size(tmp.dat,4);
        nVox = prod(size(tmp.dat,1:3));

        %% --- Preallocate y ---
        if nVolumes > 1
            assert(size(func,1)==1,'4D data only supported for single file');
            y = zeros(nVolumes,nVox,precision);
        else
            y = zeros(size(func,1),nVox,precision);
        end

        %% --- Load in batches ---
        nFiles = size(func,1);
        for i = 1:nFiles
            fprintf('.');
            if i > 1, tmp = nifti(func(i,:)); end

            if nVolumes > 1
                parfor j = 1:nVolumes
                    if strcmp(precision,'logical')
                        y(j,:) = reshape(tmp.dat(:,:,:,j),[1,nVox])~=0;
                    else
                        y(j,:) = single(reshape(tmp.dat(:,:,:,j),[1,nVox]));
                    end
                end
                break
            else
                if strcmp(precision,'logical')
                    y(i,:) = tmp.dat(:)~=0;
                else
                    y(i,:) = single(tmp.dat(:));
                end
            end
        end


        %% --- Load coordinates and affine ---
        Vfunc = spm_vol(func(1,:));
        [~,XYZmm] = spm_read_vols(Vfunc);
        XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
        dim = Vfunc.dim;
        mat = Vfunc.mat;
    end
    fprintf(' done.\n');
end