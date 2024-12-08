function afxSecondlevelPiniROIs(firstlevelInfo)

fprintf('afxSecondlevelPiniROIs ...\n');
    % defaults
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end

    % load firstlevel information
    info = load(firstlevelInfo);
    info.subjects = info.subjects(~[info.subjects.exclude]); % get rid of excluded subjects

    % define space
    [~,XYZmm,dim,mat] = afxLoadFunc(info.subjects(1).conditions(1).func);
    
    % get all conditions
    tmpCond = [info.subjects.conditions];
    conditions = unique({tmpCond.name});    
    
    % 2nd level directory and output files
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');
    dirSecondlevel = fullfile('results',info.firstlevelDir,'secondlevel','PiniROIs');
    if ~exist(dirSecondlevel,'dir'), mkdir(fullfile(dirSecondlevel,'pc1')); end
    
    % number of subjects and names
    subsAll = {info.subjects.name}';

    % to do: do a batch of ROIs at once for better performance (but not all
    % ROIs since the RAM might not be enough
    for iRoi = 1:length(info.rois)
        fprintf('   ROI %i ...',iRoi);
        % load all data
        mergedMatrix = [];
        for iSub = 1:length(subsAll)
            for iCond = 1:length(conditions)
                load(fullfile(dirFirstlevel,['cond_',conditions{iCond}],'meanWithinRoiConn',[subsAll{iSub},'.mat']))
                allDat(iSub,iCond).cc = yRoiCC(iRoi).cc;
            end
            ind = yRoiCC(iRoi).ind;
            name = yRoiCC(iRoi).name;
            % average across conditions (transposed)
            mergedMatrix(iSub,:) = mean([allDat(iSub,:).cc],2);
        end

        % remove outlier
        globalConn = mean(mergedMatrix,2);
        m = mean(globalConn);
        s = std(globalConn);
        idx = globalConn > m+3*s;
        mergedMatrix(idx,:) = [];

        % check for faulty connectivity estimates
        idx = any(isinf(mergedMatrix),2);
        if any(idx)
            mergedMatrix(idx,:) = [];
            warning([num2str(sum(idx)) ' subjects needed to be excluded due to nonreal connectivity estimates']);
        end

        % PCA
        [pCoeff] = pca(mergedMatrix);

        % project back PC1
        pc1Coeff = pCoeff(:,1);
        dat = nan(1,size(XYZmm,2));
        dat(ind) = pc1Coeff;
        % save pc1
        fname = fullfile(dirSecondlevel,'pc1',[name,'.nii']);
        afxVolumeWrite(fname,dat,dim,'int16',mat,'',true);
        % save mask (> 20th percentile)
        fname = fullfile(dirSecondlevel,[name,'.nii']);
        afxVolumeWrite(fname,dat>prctile(pc1Coeff,80),dim,'uint8',mat,'',true);
        fprintf(' done\n');
    end
    
    fprintf('afxSecondlevelPiniROIs ... done\n');
end