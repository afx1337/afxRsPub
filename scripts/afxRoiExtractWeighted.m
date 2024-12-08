function [y,names] = afxRoiExtractWeighted(rois,mask,images)
    if nargin < 1
        rois = spm_select([1 Inf],'image','Select ROIs');
    end
    for i = 1:size(rois,1)
        [~,name]  = fileparts(rois(i,:));
        tmp(i).name   = name;
        tmp(i).file   = rois(i,:);
    end
    rois = tmp;
        
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
        dat = afxVolumeResample(curRoi.file,XYZmm,1);
        
        y(iRoi).name = curRoi.name;
        y(iRoi).z = atanh(corr(X(:,mask)',dat(mask)));
    end
    names = {Vi.fname}';
    
    out = [['-';names] [{y.name};num2cell([y.z])]];
    [filename, pathname] = uiputfile('*.xlsx', 'Save as ...');
    if ischar(filename)
        destFile = fullfile(pathname,filename);
        xlswrite(destFile,out);
    end
end