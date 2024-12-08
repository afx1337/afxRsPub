function afxNetworkExtract(fNetwork,fImages)
    if nargin < 1
        fNetwork = spm_select([1 1],'image','Select network');
    end
    if nargin < 2
        fImages = spm_select([1 Inf],'image','Select images');
    end
    pm  = spm_input(['premultiply by'],1,'r',1,1);
    u  = spm_input(['Hight threshold (0 for none)'],'+1','r',0,1);
    k  = spm_input(['Extend threshold (0 for none)'],'+1','r',0,1);

    outVals = [];
    %sizes = [];
    outRoiNames = {};
  
    % load images
    Vi = spm_vol(char(fImages));
    n = length(Vi);
    parfor iRoi = 1:n
        tmp = spm_read_vols(Vi(iRoi));
        val(iRoi,:) = tmp(:);
    end
    outimageNames = {Vi.fname}';
    Vi = Vi(1);
    [~,XYZmm] = spm_read_vols(Vi);
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    
    % load network
    dat = afxVolumeResample(fNetwork,XYZmm,0)*pm;
    ind = find(dat > u);
    [X,Y,Z] = ind2sub(Vi.dim,ind);
    XYZ = [X Y Z]';
    Z = dat(ind);
    val = val(:,ind);
    clear dat X Y;
    
    % extract data per cluster
    A = spm_clusters(XYZ);
    for iRoi=1:numel(unique(A))
        if sum(A==iRoi)>=k
            nXYZ  = XYZ(:,A==iRoi);
            mXYZ  = mean(Vi.mat * [nXYZ; ones(1,size(nXYZ,2))],2); mXYZ = mXYZ(1:3)';

            fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
            spm_orthviews('Image', spm_vol(fullfile(spm('dir'),'canonical','single_subj_T1.nii')), [0.0 0.22 1 .8]);
            spm_orthviews('addcolouredblobs',1,XYZ,Z,Vi(1).mat,[.5 0 0]);
            spm_orthviews('addcolouredblobs',1,nXYZ,Z(A==iRoi),Vi(1).mat,[1 1 0]);
            spm_orthviews('reposition',mXYZ)

            ROIname = spm_input('VOI name [0 skips]','+1','s',[spm_str_manip(fNetwork,'rt') '_' num2str(iRoi)]);

            if ~strcmp(ROIname,'0')
                outRoiNames{end+1} = ROIname;
                outVals(1:size(val,1),end+1) = mean(val(:,A==iRoi),2);
                %sizes(1,end+1) = sum(~isnan(val(1,A==iRoi)));
            end
          end

    end
    
    % image names
    for i = 1:size(outimageNames,1)
        [p,f,e] = fileparts(outimageNames{i});
        n1{i} = p;
        n2{i} = strcat(f,e);
    end
    
    % save as xlsx
    %out = [ 'imagePath' 'imageName' outRoiNames ; '-' 'size' num2cell(sizes) ; n1' n2'  num2cell(outVals)];
    out = [ 'imagePath' 'imageName' outRoiNames ; n1' n2'  num2cell(outVals)];
    [filename, pathname] = uiputfile('*.xlsx', 'Save as ...');
    if ischar(filename)
        destFile = fullfile(pathname,filename);
        xlswrite(destFile,out);
    end
end