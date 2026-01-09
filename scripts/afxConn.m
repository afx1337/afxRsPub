function rois = afxConn(func,masks,rpFile,denoisingOptions,rois,firstlevelDir,conditionName,subjectName,analyses)

    % output directory
    outDir = fullfile('results',firstlevelDir,'firstlevel',['cond_' conditionName]);

    % load functional data (smoothed/unsmoothed)
    if denoisingOptions.unsmoothed
        [y,XYZmm,dim,mat] = afxLoadFunc(func.func2);
    else
        [y,XYZmm,dim,mat] = afxLoadFunc(func.func);
    end
    
    % load masks (brainmask is hardcoded)
    [brainMask,subjectMasks] = afxLoadMasks(masks,fullfile('masks','brainmask.nii'),XYZmm);
    
    % denoising (confound removal, filtering, scrubbing)
    y = afxDenoising(y,subjectMasks,brainMask,rpFile,denoisingOptions);
    y(:,~brainMask) = NaN;
    
    % roi timeseries extraction
    if ~isempty(rois) && any(contains(analyses,{'fc_wholebrain','fc_network'}))
        [yRoi, lesioned] = afxRoiTimeseries(rois,y,subjectMasks,brainMask,XYZmm,denoisingOptions);
    end
    
    % smoothing
    if denoisingOptions.unsmoothed
        y = afxSmoooth(y,denoisingOptions.sFWHM,dim,mat);
    end
    
    % calculation of functional connectivity
    for i = 1:length(analyses)
        switch analyses{i}
            case 'quality'
                afxQuality(y,subjectMasks,dim,mat,outDir,subjectName);
            case 'fc_wholebrain'
                afxConnWholeBrain(y,yRoi,brainMask,rois,dim,mat,outDir,subjectName);
            case 'fc_network'
                afxConnNetwork(yRoi,lesioned,rois,outDir,subjectName);
            case 'vmhc'
                afxConnVMHC(y,brainMask,dim,mat,outDir,subjectName);
            case 'hemo_lag'
                afxConnHemodynamicLag(y,brainMask,subjectMasks,rpFile,denoisingOptions,dim,mat,outDir,subjectName);
            case 'fALFF'
                afxfALFF(y,brainMask,subjectMasks,denoisingOptions,dim,mat,outDir,subjectName);
            case 'ALFF'
                afxALFF(y,brainMask,subjectMasks,denoisingOptions,dim,mat,outDir,subjectName);                
            case 'piniroi'
                afxConnPiniRoi(y,brainMask,subjectMasks,rois,XYZmm,denoisingOptions,outDir,subjectName);
            otherwise
                error(['Unknown analysis: ' analyses{i}]);
        end
    end
end