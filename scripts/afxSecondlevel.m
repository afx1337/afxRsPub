function afxSecondlevel(firstlevelInfo,covariateFactors,groupFactor,conditions,mode)
    spm('Defaults','fmri');
    spm_jobman('initcfg');

    % defaults
    if nargin < 5, mode = 'roi'; end
    if nargin < 4, conditions = []; end
    if nargin < 3, groupFactor = []; end
    if nargin < 2, covariateFactors = {}; end
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end

    % load firstlevel information
    info = load(firstlevelInfo);
    
    % create pseudo roi if non-roi based approach
    if ~strcmp(mode,'roi')
        info.rois = struct;
        info.rois(1).name = mode;
    end
        
    % make 2nd level directory
    if isempty(groupFactor), tmpG = 'all'; else tmpG = groupFactor; end
    if isempty(covariateFactors), tmpCov = 'none'; else tmpCov = [covariateFactors{:}]; end
    if isempty(conditions), tmpCond = 'all'; else tmpCond = [conditions{:}]; end
    name = sprintf('%s_group-%s_cond-%s_cov-%s',mode,tmpG,tmpCond,tmpCov);
    dirSecondlevel = fullfile('results',info.firstlevelDir,'secondlevel','flexibleFactorial',name);
    fSecondlevel = fullfile(dirSecondlevel,'SPM.mat');
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');
    mkdir(dirSecondlevel);
    
    % get all conditions
    if isempty(conditions)
       tmpCond = [info.subjects(~[info.subjects.exclude]).conditions];
       conditions = unique({tmpCond.name});
    end
    
    % define matlabbatch structure for 2nd level
    matlabbatch{1}.spm.stats.factorial_design.dir = { dirSecondlevel };
    
    % define factors for factorial design
    iFactor = 1;
    factors = struct('subject',1,'condition',[],'group',[],'roi',[]);
    % subject factor
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).name = 'subject';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).ancova = 0;
    iFactor = iFactor + 1;
    % factor group (between-subject) factor
    if ~isempty(groupFactor)
        factors.group = iFactor;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).name = 'group';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).ancova = 0;
        iFactor = iFactor + 1;
        tmp = [info.subjects(~[info.subjects.exclude]).info];
        if ischar(tmp(1).(groupFactor))
            groups = unique({tmp.(groupFactor)});
        else
            groups = unique([tmp.(groupFactor)]);
        end
    end
    % factor condition (within-subject) factor
    if length(conditions) > 1
        factors.condition = iFactor;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).name = 'condition';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).dept = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).ancova = 0;
        iFactor = iFactor + 1;
    end
    % factor roi (within-subject) factor
    factors.roi = iFactor;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).name = 'roi';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(iFactor).ancova = 0;

    % initialise covariates
    cov = struct([]);
    for iCov = 1:length(covariateFactors)
        cov(iCov).name = covariateFactors{iCov};
        cov(iCov).dat = [];
    end
    
    % add scans, conditions and covariates
    subInd = find(~[info.subjects.exclude]);
    iSubjectReal = 1;
    subjectsInclude = false(1,length(subInd));
    for iSubject = 1:length(subInd)
        subject = info.subjects(subInd(iSubject));
        if ~isempty(factors.group)
            try
                group(iSubject) = find(strcmp(groups,subject.info.(groupFactor)));
            catch
                group(iSubject) = find(groups==subject.info.(groupFactor));
            end
        end
        iImg = 1;
        conds = [];
        for iCond = 1:length(subject.conditions)
            for iRoi = 1:length(info.rois)
                % add imaage
                if strcmp(mode,'fALFF') || strcmp(mode,'ALFF')
                    imgFirstlevel =  fullfile(pwd,dirFirstlevel,['cond_' subject.conditions(iCond).name],mode,[subject.name '.nii']);
                else
                    imgFirstlevel =  fullfile(pwd,dirFirstlevel,['cond_' subject.conditions(iCond).name],['roi_' info.rois(iRoi).name],[subject.name '.nii']);
                end
                if exist(imgFirstlevel,'file') && any(strcmp(conditions,subject.conditions(iCond).name))
                    subjectsInclude(iSubject) = true;
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(iSubjectReal).scans{iImg,1} = imgFirstlevel;
                    % get conditions of current image
                    newCondRow = nan(1,iFactor);
                    if ~isempty(factors.group)
                        newCondRow(1,factors.group) = group(iSubject);
                    end
                    if ~isempty(factors.condition)
                        newCondRow(1,factors.condition) = find(strcmp(conditions,subject.conditions(iCond).name));
                    end
                    if ~isempty(factors.roi)
                        newCondRow(1,factors.roi) = iRoi;
                    end
                    conds = [conds; newCondRow(1,2:end)];
                    % get covariates for current image
                    for iCov = 1:length(covariateFactors)
                        if strcmp(covariateFactors{iCov},'movement')
                            cov(iCov).dat(end+1) = subject.conditions(iCond).movement.meanFD;
                        else
                            cov(iCov).dat(end+1) = subject.info.(covariateFactors{iCov});
                        end
                    end
                    iImg = iImg+1;
                end
            end
        end
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(iSubjectReal).conds = conds;
        iSubjectReal = iSubjectReal + 1;
    end
    
    % which factors (and interactions) shall be included in the model
    % main factors
    % include subject factor only for small sample sizes
    iFacInt = 1;
    if sum(~[info.subjects.exclude]) < 70 && isempty(factors.group)
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.fmain.fnum = 1;
        nSub = nnz(subjectsInclude);
        iFacInt = iFacInt + 1;
    else
        nSub = 0;
    end
    % group factor
    if ~isempty(factors.group)
        group = group(subjectsInclude);
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.fmain.fnum = factors.group;
        for iGroup = 1:length(groups)
            weightsGroup(iGroup) = nnz(group == iGroup)/length(group);
            if nSub ~= 0
                subsGroup(iGroup,:) = (group == iGroup);
            end
        end
        if nSub == 0, subsGroup = zeros(length(groups),0); end
        iFacInt = iFacInt + 1;
    end
    % condition factor
    if ~isempty(factors.condition)
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.fmain.fnum = factors.condition;
        iFacInt = iFacInt + 1;
    end
    % roi factor
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.fmain.fnum = factors.roi;
    iFacInt = iFacInt + 1;
    
    % interactions (only two way interactions ...)
    if isempty(factors.group) && ~isempty(factors.condition)
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.inter.fnums = [factors.condition factors.roi];
        iFacInt = iFacInt + 1;
    end
    if ~isempty(factors.group) && isempty(factors.condition)
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{iFacInt}.inter.fnums = [factors.group factors.roi];
        iFacInt = iFacInt + 1;
    end
    
    
    % fill covariates to matlabbatch
    if isempty(covariateFactors)
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    else
        for iCov = 1:length(covariateFactors)
            matlabbatch{1}.spm.stats.factorial_design.cov(iCov).c = cov(iCov).dat;
            matlabbatch{1}.spm.stats.factorial_design.cov(iCov).cname = cov(iCov).name;
            matlabbatch{1}.spm.stats.factorial_design.cov(iCov).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(iCov).iCC = 1;
        end
    end

    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = { fullfile(pwd,'masks','brainmask.nii') };
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %%%% CONTRASTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nRois = length(info.rois);
    nCov = length(cov);
    nCond = length(conditions);
    
    % no group and no condition factor => only roi main effects
    if isempty(factors.group) && isempty(factors.condition)
        tmpRois = eye(nRois);
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        iCon = 1;
        for iRoi = 1:nRois
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['ROI: ' info.rois(iRoi).name '+'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ tmpRois(iRoi,:) zeros(1,nCov) ones(1,nSub)/nSub];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['ROI: ' info.rois(iRoi).name '-'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = -[ tmpRois(iRoi,:) zeros(1,nCov) ones(1,nSub)/nSub];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
        end
        matlabbatch{3}.spm.stats.con.delete = 1;
    end

    % no group factor => 1) roi main effect, 2) condition x roi interaction,
    % 3) condition differences x roi interaction
    if isempty(factors.group) && ~isempty(factors.condition)
        tmpRois = eye(nRois);
        tmpCond = eye(nCond);
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        % roi main effects
        iCon = 1;
        for iRoi = 1:nRois
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['main effect ROI ' info.rois(iRoi).name '+'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ ones(1,nCond)/nCond tmpRois(iRoi,:) repmat(tmpRois(iRoi,:),1,nCond)/nCond zeros(1,nCov) ones(1,nSub)/nSub];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['main effect ROI ' info.rois(iRoi).name '-'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = -[ ones(1,nCond)/nCond tmpRois(iRoi,:) repmat(tmpRois(iRoi,:),1,nCond)/nCond zeros(1,nCov) ones(1,nSub)/nSub];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
        end
        % condition x roi interaction
        for iCond = 1:nCond
            for iRoi = 1:nRois
                interaction = tmpRois(iRoi,:)'*tmpCond(iCond,:);
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Cond ' conditions{iCond} ' x ROI ' info.rois(iRoi).name '+'];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ tmpCond(iCond,:) tmpRois(iRoi,:) interaction(:)' zeros(1,nCov) ones(1,nSub)/nSub];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
                interaction = tmpRois(iRoi,:)'*tmpCond(iCond,:);
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Cond ' conditions{iCond} ' x ROI ' info.rois(iRoi).name '-'];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = -[ tmpCond(iCond,:) tmpRois(iRoi,:) interaction(:)' zeros(1,nCov) ones(1,nSub)/nSub];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
            end
        end
        % condition differences x roi interaction
        [cond1,cond2] = find(1-eye(nCond));
        for iDiff = 1:length(cond1)
            tmpCond = zeros(1,nCond);
            tmpCond(cond1(iDiff)) = 1;
            tmpCond(cond2(iDiff)) = -1;
            for iRoi = 1:nRois
                interaction = tmpRois(iRoi,:)'*tmpCond;
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Cond ' conditions{cond1(iDiff)} ' > ' conditions{cond2(iDiff)} ' x ROI ' info.rois(iRoi).name];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ tmpCond zeros(1,nRois) interaction(:)' zeros(1,nCov) ];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
            end
        end
        matlabbatch{3}.spm.stats.con.delete = 1;
    end

    % no condition factor => 1) roi main effect, 2) group x roi interaction,
    % 3) group differences x roi interaction
    if ~isempty(factors.group) && isempty(factors.condition)
        nGroups = length(groups);
        tmpRois = eye(nRois);
        tmpGroup = eye(nGroups);
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        % roi main effects
        iCon = 1;
        for iRoi = 1:nRois
            tmp = tmpRois(iRoi,:)'*weightsGroup;
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['main effect ROI ' info.rois(iRoi).name '+'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ weightsGroup tmpRois(iRoi,:) tmp(:)' zeros(1,nCov) ones(1,nSub)/nSub ];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['main effect ROI ' info.rois(iRoi).name '-'];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = -[ weightsGroup tmpRois(iRoi,:) tmp(:)' zeros(1,nCov) ones(1,nSub)/nSub ];
            matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
            iCon = iCon + 1;
        end
        % group x roi interaction
        % group labels to strings
        if isnumeric(groups)
            for i = 1:length(groups)
                groupsNew{i} = num2str(groups(i));
            end
            groups = groupsNew;
        end
        for iGroup = 1:nGroups
            for iRoi = 1:nRois
                interaction = tmpRois(iRoi,:)'*tmpGroup(iGroup,:);
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Group ' groups{iGroup} ' x ROI ' info.rois(iRoi).name '+'];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ tmpGroup(iGroup,:) tmpRois(iRoi,:) interaction(:)' zeros(1,nCov) subsGroup(iGroup,:)/nnz(subsGroup(iGroup,:)) ];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Group ' groups{iGroup} ' x ROI ' info.rois(iRoi).name '-'];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = -[ tmpGroup(iGroup,:) tmpRois(iRoi,:) interaction(:)' zeros(1,nCov) subsGroup(iGroup,:)/nnz(subsGroup(iGroup,:)) ];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
            end
        end
        % group differences x roi interaction
        [gr1,gr2] = find(1-eye(nGroups));
        for iDiff = 1:length(gr1)
            tmpGroup = zeros(1,nGroups);
            tmpGroup(gr1(iDiff)) = 1;
            tmpGroup(gr2(iDiff)) = -1;
            for iRoi = 1:nRois
                interaction = tmpRois(iRoi,:)'*tmpGroup;
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.name    = ['Group ' groups{gr1(iDiff)} ' > ' groups{gr2(iDiff)} ' x ROI ' info.rois(iRoi).name];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.weights = [ tmpGroup zeros(1,nRois) interaction(:)' zeros(1,nCov) ];
                matlabbatch{3}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
                iCon = iCon + 1;
            end
        end
        matlabbatch{3}.spm.stats.con.delete = 1;
    end
    % save and run 2nd level batch
    save(fullfile(dirSecondlevel,'2nd_level_batch.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
end