function [y,names] = afxRoiExtractMasked(rois,mask,images)
    if nargin < 1
        tmp = spm_select([1 Inf],'image','Select ROIs');
        for i = 1:size(tmp,1)
            [~,name]  = fileparts(tmp(i,:));
            rois(i).type   = 'image';
            rois(i).name   = name;
            rois(i).file   = tmp(i,:);
        end
    end
    
    if nargin < 2
         mask = spm_select([1 1],'image','Select mask');
    end
        
    if nargin < 3
        images = spm_select([1 Inf],'image','Select images');
    end
    
    Vi = spm_vol(char(images));
    n = length(Vi);
    parfor i = 1:n
        tmp = spm_read_vols(Vi(i));
        X(i,:) = tmp(:);
    end
    [~,XYZmm] = spm_read_vols(Vi(1));
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    mask = afxVolumeResample(mask,XYZmm,0) > .5;

    y = struct([]);
    for iRoi = 1:length(rois)
        % get ROI voxel indices
        curRoi = rois(iRoi);
        switch curRoi.type
            case 'sphere'
                xyz = curRoi.coords;
                D = sqrt((XYZmm(1,:)-xyz(1)).^2+(XYZmm(2,:)-xyz(2)).^2+(XYZmm(3,:)-xyz(3)).^2);
                ind = (D < curRoi.radius)';
                clear D;
            case 'image'
                dat = afxVolumeResample(curRoi.file,XYZmm,0);
                ind = dat > 0.5;
                clear dat;
            otherwise
                error('unknown roi type: %s',curRoi.type);
        end
        
        y(iRoi).name = curRoi.name;
        y(iRoi).y = X(:,ind & mask);
        y(iRoi).sizeTotal = sum(ind);
        y(iRoi).sizeMask = sum(ind & mask);
        y(iRoi).sizeRatio = sum(ind & mask)/sum(ind);
        y(iRoi).mean = mean(X(:,ind & mask),2);
    end
    names = {Vi.fname}';
end