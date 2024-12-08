function [allMean,allMeanZ,allT] = afxLnmSecondlevel(firstlevelInfo,groups,deleteFirstlevel)
    if nargin < 3, deleteFirstlevel = false; end
    if nargin < 2, groups = []; end
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end

    % load firstlevel information
    info = load(firstlevelInfo);
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');
    dirSecondlevel = fullfile('results',info.firstlevelDir,'secondlevel','lnsm');
    
    spm_jobman('initcfg'); % afxNlmMean uses jobmanager
    
    allT = cell(length(info.rois),1);
    allMean = cell(length(info.rois),1);
    allMeanZ = cell(length(info.rois),1);
    
    % mean ROI connectivity
    for iRoi = 1:length(info.rois)
        fil = {};
        for iSubject =  find(~[info.subjects.exclude])
            for iCond = 1:length(info.subjects(iSubject).conditions)
                fil{end+1,1} = fullfile(dirFirstlevel,['cond_' info.subjects(iSubject).conditions(iCond).name],['roi_' info.rois(iRoi).name],[info.subjects(iSubject).name '.nii']);
            end
        end
        
        fMean  = fullfile(dirSecondlevel,'mean',['mean_roi_' info.rois(iRoi).name '.nii']);
        %fMeanZ = fullfile(dirSecondlevel,'meanZ',['meanZ_roi_' info.rois(iRoi).name '.nii']);
        fT  = fullfile(dirSecondlevel,'ttest',['tmap_roi_' info.rois(iRoi).name '.nii']);
        
        afxLnmMean(fil,fMean);
        %afxLnmZ(fMean,fMeanZ)
        afxLnmTTest(fil,fT,iCond);
        
        allT{iRoi} = fT;
        allMean{iRoi} = fMean;
        %allMeanZ{iRoi} = fMeanZ;

        % delete firstlevel images
        if deleteFirstlevel
            for iFile = 1:length(fil)
                delete(fil{iFile});
            end
        end
    end

    % delete firstlevel directory
    if deleteFirstlevel
        rmdir(dirFirstlevel);
    end
    
    % Threshold .00005, see Boes et al., 2015
    alpha = .00005;
    df = length(find(~[info.subjects.exclude]))-1;
    tCrit = tinv(1-alpha,df);
    if ~isempty(groups)
        for i = 1:length(groups)
            afxLnmSum(allT(groups{i}),fullfile(dirSecondlevel,'ttest',['ttest_pos_sum_group' num2str(i) '.nii']),tCrit);
            afxLnmSum(allT(groups{i}),fullfile(dirSecondlevel,'ttest',['ttest_neg_sum_group' num2str(i) '.nii']),-tCrit);
        end
        afxLnmTTest2(allMean(groups{2}),allMean(groups{1}),.001,fullfile(dirSecondlevel,'mean','ttest_group2_vs_group1.nii'));
        %afxLnmTTest2(allMeanZ(groups{2}),allMeanZ(groups{1}),.001,fullfile(dirSecondlevel,'meanZ','ttest_group2_vs_group1.nii'));
    end
    fTCrit = fullfile(dirSecondlevel,'ttest','info.txt'); 
    fid = fopen(fTCrit,'wt');
    fprintf(fid, 'df = %i\nwith alpha = %f and one-tailed testing, tCrit = %f\n',df,alpha,tCrit);
    fclose(fid);
end
