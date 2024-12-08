function afxMotionStats()
    subjectsFile = spm_select(1,'.mat','Select file ...',{},'data','^subjects*.mat');
    [pth,fname,~] = fileparts(subjectsFile);
    statsFile = fullfile(pth,strcat(fname,'-demographics.xlsx'));
    load(subjectsFile);

    out = cell(0);
    out(end+1,:) = { 'Name' 'Condition' 'FD > .5' 'Mean(FD)' 'Exclusion' 'Reason' 'Age' 'Sex (1=male)' 'Handedness (1=right)' };

    for iSub = 1:length(subjects)
        outName = subjects(iSub).name;
        for iCond = 1:length(subjects(iSub).conditions)
            curCond = subjects(iSub).conditions(iCond);
            outCond = curCond.name;
            outFd5 = nnz(curCond.movement.FD > .5);
            outMeanFd = mean(curCond.movement.FD);
            outExclusion = subjects(iSub).exclude;
            outReason = '';
            if subjects(iSub).exclude, outReason = strjoin(subjects(iSub).excludeReason,'; '); end
            if isfield(subjects(iSub).info,'age'), outAge = subjects(iSub).info.age; else, outAge = ''; end
            if isfield(subjects(iSub).info,'sex'), outSex = subjects(iSub).info.sex; else, outSex = ''; end
            if isfield(subjects(iSub).info,'handedness'), outHandedness = subjects(iSub).info.handedness; else, outHandedness = ''; end
            out(end+1,:) = {outName outCond outFd5 outMeanFd outExclusion outReason outAge outSex outHandedness };
        end
    end
    xlswrite(statsFile,out,1);
    xlswrite(statsFile,out(logical([1 [out{2:end,5}] == 0]),:),2);
end