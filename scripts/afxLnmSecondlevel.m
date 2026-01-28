function [allMean,allT] = afxLnmSecondlevel(firstlevelInfo,groups,deleteFirstlevel)
    if nargin < 3, deleteFirstlevel = false; end
    if nargin < 2, groups = []; end
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end

    % load firstlevel information
    info = load(firstlevelInfo);
    df = length(find(~[info.subjects.exclude]))-1;
    nCond = length(info.subjects(1).conditions); % all subjects should have same number of conditions
    tCrit = tinv(1-0.001,df);
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');
    dirSecondlevel = fullfile('results',info.firstlevelDir,'secondlevel','lnsm');
    mkdir(fullfile(dirSecondlevel,'mean'));
    mkdir(fullfile(dirSecondlevel,'ttest'));
    mkdir(fullfile(dirSecondlevel,'ttest_img'));
    
    allT = cell(length(info.rois),1);
    allMean = cell(length(info.rois),1);
    
    % for each ROI
    for iRoi = 1:length(info.rois)
        fil = {};
        for iSubject =  find(~[info.subjects.exclude])
            for iCond = 1:nCond
                fil{end+1,1} = fullfile(dirFirstlevel,['cond_' info.subjects(iSubject).conditions(iCond).name],['roi_' info.rois(iRoi).name],[info.subjects(iSubject).name '.nii']);
            end
        end
        
        fMean  = fullfile(dirSecondlevel,'mean',['mean_roi_' info.rois(iRoi).name '.nii']);
        fT = fullfile(dirSecondlevel,'ttest',['tmap_roi_' info.rois(iRoi).name '.nii']);
        fTPng = fullfile(dirSecondlevel,'ttest_img',['tmap_roi_' info.rois(iRoi).name '.png']);
        
        afxLnmSecondlevelCalc(fil, fMean, fT, nCond);
        
        % create sections to allow reviewing individual lesion networks
        if ~exist(fTPng,'file')
            f = afxPlotSections(fT,[40 20 0 -30],{},[tCrit Inf],[tCrit Inf],1,3);
            set(f,'Position',[100 100 1500/2 360/2]);
            print(f,fTPng,'-dpng','-r150');
            close(f);
        end
        
        allT{iRoi} = fT;
        allMean{iRoi} = fMean;
    end

    % delete firstlevel directory
    % note, that this deletes *all* firstlevel analyses
    if deleteFirstlevel
        rmdir(dirFirstlevel,'s');
    end
    
    % Boes threshold .00005, (Boes et al., 2015)
    alphaBoes = .00005;
    tCritBoes = tinv(1-alphaBoes,df);
    
    % save tCrit to txt file
    fTCrit = fullfile(dirSecondlevel,'ttest','info.txt'); 
    fid = fopen(fTCrit,'wt');
    fprintf(fid, 'df = %i\nwith alpha = %f and one-tailed testing, tCrit = %f\n',df,alphaBoes,tCritBoes);
    fprintf(fid, 'alpha = 0.001 corresponds to tCrit = %f\n',tCrit);
    fclose(fid);
    
    % legacy secondlevel code used for comparibility with Boes et al., 2015
    %if ~isempty(groups)
    %    spm_jobman('initcfg');
    %    for i = 1:length(groups)
    %        afxLnmSum(allT(groups{i}),fullfile(dirSecondlevel,'ttest',['ttest_pos_sum_group' num2str(i) '.nii']),tCritBoes);
    %        afxLnmSum(allT(groups{i}),fullfile(dirSecondlevel,'ttest',['ttest_neg_sum_group' num2str(i) '.nii']),-tCritBoes);
    %    end
    %    afxLnmTTest2(allMean(groups{2}),allMean(groups{1}),.001,fullfile(dirSecondlevel,'mean','ttest_group2_vs_group1.nii'));
    %end
end
