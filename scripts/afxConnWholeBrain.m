function afxConnWholeBrain(y,yRoi,brainMask,rois,dim,mat,outDir,subjectName)
    % calculation of functional connectivity
        
    fprintf('   Calculating whole brain functional connectivity ...\n')
    
    % prepare whole brain time series
    nVox = size(y,2);
    y = y(:,brainMask);           % keep only data within the brain mask
    y = y - mean(y,1);            % mean centering
    y = y ./ vecnorm(y);          % norm all data
    
    % prepare ROI time series
    yRoi = yRoi - mean(yRoi,1);   % mean centering
    yRoi = yRoi ./ vecnorm(yRoi); % norm roi timeseries
    
    % blockwise (speedup BS 25 -> 7-8x)
    blockSize = max(1,round(250000/nVox*25)); % constant memory usage
    nRoi = size(yRoi,2);
    roiIdx = 1:blockSize:nRoi;  % start indices of each block    
    for b = 1:length(roiIdx)
        idxStart = roiIdx(b);
        idxEnd   = min(idxStart + blockSize - 1, nRoi);
        fprintf('      Processing ROIs %d-%d ... ', idxStart, idxEnd);
        idxBlock = idxStart:idxEnd;

        % compute fc
        Zblock = yRoi(:,idxBlock)' * y;  % size: blockSize × nVox

        fprintf('saving to disk ')

     
        for i = 1:length(idxBlock)
            iRoi = idxBlock(i);

            % fisher transformation
            z = min(Zblock(i,:), 1-1e-15); % to prevent this: atanh(1) = Inf
            z = atanh(z);

            % voxels outside brain mask are NaN
            z2 = nan(1,nVox,'like',y);
            z2(brainMask) = z;
            clear z;

            % create output dir if necassary
            imgFname = fullfile(outDir,['roi_' rois(iRoi).name],strcat(subjectName,'.nii'));
            [pth,~,~] = fileparts(imgFname);
            if ~exist(pth,'dir'), mkdir(pth); end

            % save fucnctional connectivity map
            afxVolumeWrite(imgFname,z2,dim,'int16',mat,'fisher transformed corrcoeff');
            fprintf('.')
        end
        fprintf(' done\n')
    end
end
