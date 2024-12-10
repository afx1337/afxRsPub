function afxConnWholeBrain(y,yRoi,brainMask,rois,dim,mat,outDir,subjectName)
    % calculation of functional connectivity
    fprintf('   Calculating whole brain functional connectivity ...\n')
    
    % prepare whole brain time series
    % keep only data within the brain mask (saves memory in the parfor-loop)
    nVox = size(y,2);
    y = y(:,brainMask);
    % mean centering
    y = y - repmat(mean(y),size(y,1),1);
    % norm all data
    parfor i = 1:size(y,2)
        y(:,i) = y(:,i) ./ norm(y(:,i));
    end
    
    % prepare ROI time series
    % mean centering
    yRoi = yRoi - repmat(mean(yRoi),size(yRoi,1),1);
    % norm roi timeseries
    parfor iRoi = 1:size(yRoi,2)
        yRoi(:,iRoi) = yRoi(:,iRoi) ./ norm(yRoi(:,iRoi));
    end
    
    for iRoi = 1:size(yRoi,2)
        fprintf('      > Roi %s ... ',rois(iRoi).name);
        % compute fc
        y2 = yRoi(:,iRoi)';
        z = nan(1,size(y,2));
        parfor i = 1:size(y,2)
            % a = a - mean(a); b = b - mean(b)
            % rho = a*b / ( norm(a)*norm(b) )
            % note, that mean centering was done before for roi timeseries and
            % is not necessary for voxel timeseries due to highpass filtering
            % roi timeseries is already normed
            z(i) =  y2 * y(:,i);
        end
        % fisher transformation
        z(z > 1-1e-15) = 1-1e-15; % to prevent this: atanh(1) = Inf
        z = atanh(z);

        % voxels outside brain mask are NaN
        z2 = nan(1,nVox);
        z2(brainMask) = z;
        clear z;
        
        % create output dir if necassary
        imgFname = fullfile(outDir,['roi_' rois(iRoi).name],strcat(subjectName,'.nii'));
        [pth,~,~] = fileparts(imgFname);
        if ~exist(pth,'dir'), mkdir(pth); end
        % save fucnctional connectivity map
        afxVolumeWrite(imgFname,z2,dim,'int16',mat,'fisher transformed corrcoeff');
        fprintf('done\n')
    end
end
