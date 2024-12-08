function afxMeanFirstlevels(firstlevelInfo)
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end
    
    info = load(firstlevelInfo);
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');

    % mean ROI connectivity
    for iRoi = 1:length(info.rois)
        for iSubject =  find(~[info.subjects.exclude])
            fil = {};
            for iCond = 1:length(info.subjects(iSubject).conditions)
                fil{end+1,1} = fullfile(dirFirstlevel,['cond_' info.subjects(iSubject).conditions(iCond).name],['roi_' info.rois(iRoi).name],[info.subjects(iSubject).name '.nii']);
            end
            
            fMean  = fullfile(dirFirstlevel,'mean',['roi_' info.rois(iRoi).name],['mean_' info.subjects(iSubject).name '.nii']);
            afxLnmMean(fil,fMean);
        end
    end
end