function afxConnPiniRoi(y,brainMask,subjectMasks,rois,XYZmm,denoisingOptions,outDir,subjectName)
    % calculation of functional connectivity
    fprintf('   Calculating within ROI functional connectivity (Pini et al., 2022) ...\n')
    
    % extract time series for all rois
    [~, ~, ~, yRoiCC] = afxRoiTimeseries(rois,y,subjectMasks,brainMask,XYZmm,denoisingOptions);
   
    % create output dir if necassary
    outFname = fullfile(outDir,'meanWithinRoiConn',strcat(subjectName,'.mat'));
    [pth,~,~] = fileparts(outFname);
    if ~exist(pth,'dir'), mkdir(pth); end
    % save data
    save(outFname,'yRoiCC');
    fprintf('done\n')
end
