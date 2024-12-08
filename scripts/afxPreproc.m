function [subjects] = afxPreproc(func,struc,options,cond,sampleName,subjectName,info)
    
    if nargin < 7, info = struct([]); end

    outDir = fullfile(pwd,'data',sampleName,'subjects',subjectName);
    if exist(fullfile(outDir,'func'),'dir')
        disp(['Skipping subject >' subjectName '< (preprocessing has already been done)']);
        return;
    end
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % rename functional files which start with "w"
    for i = 1:length(func)
        for j = 1:length(func{i})
            [tmpP,tmpN,tmpE] = fileparts(func{i}{j});
            if strcmpi(tmpN(1),'w')
                newName = fullfile(tmpP,['o' tmpN tmpE]);
                movefile(func{i}{j},newName);
                func{i}{j} = newName;
            end
        end
    end
    
    % discard dummy scans
    for i = 1:length(func), func{i} = func{i}(options.dummmyScans+1:end); end
   
    % init SPM
    spm('Defaults','fmri');
    spm_jobman('initcfg');

    % if lesion is present, set lesion to zero in structural image
    % (according to post from John Ashburner on SPM mailing list)
    if isfield(options,'lesion')
        fprintf('Perform lesion masking on structural image\n');
        clear matlabbatch;
        [strucPth,strucName,strucExt] = fileparts(struc);
        strucMasked = strcat(strucName,'-lesion_masked',strucExt);
        matlabbatch{1}.spm.util.imcalc.input = {
                                                struc
                                                options.lesion
                                                };
        matlabbatch{1}.spm.util.imcalc.output = strucMasked;
        matlabbatch{1}.spm.util.imcalc.outdir = { strucPth };
        matlabbatch{1}.spm.util.imcalc.expression = '(i2<.5).*i1';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        % run batch
        spm_jobman('run',matlabbatch);
        % set struc to new masked struct image
        struc = fullfile(strucPth,strucMasked);
    end
    
    % (optionally) perform slice time correction
    if isfield(options,'sliceorder') && ~isempty(options.sliceorder)
        fprintf('Slice time correction\n');
        clear matlabbatch;
        matlabbatch{1}.spm.temporal.st.scans = func;
        matlabbatch{1}.spm.temporal.st.nslices = length(options.sliceorder);
        matlabbatch{1}.spm.temporal.st.tr = options.TR;
        matlabbatch{1}.spm.temporal.st.ta = options.TR-(options.TR/length(options.sliceorder));
        matlabbatch{1}.spm.temporal.st.so = options.sliceorder;
        matlabbatch{1}.spm.temporal.st.refslice = round((max(options.sliceorder)-min(options.sliceorder))/2 + min(options.sliceorder));
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        % save slice timing batch
        save(fullfile(outDir,'batch_slicetiming.mat'),'matlabbatch');
        % run batch
        spm_jobman('run',matlabbatch);
        % update functional filenames (a*)
        for i = 1:length(func)
            for j = 1:length(func{i})
                [tmpP,tmpN,tmpE] = fileparts(func{i}{j});
                func{i}{j} = fullfile(tmpP,['a' tmpN tmpE]);
            end
        end
    end
    
    % define matlabbatch structure for preprocessing
    clear matlabbatch;
    % files (func and struc)
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'func';
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = func;
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'struc';
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{struc}};
    % realignment (func)
    for i = 1:length(func)
        matlabbatch{3}.spm.spatial.realign.estwrite.data{i}(1) = cfg_dep(strcat('Named File Selector: func(',num2str(i),') - Files'), substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{i}));
    end
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    % coregister (struct -> func)
    matlabbatch{4}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{4}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Named File Selector: struc(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    if isfield(options,'lesion')
        matlabbatch{4}.spm.spatial.coreg.estimate.other = { options.lesion };
    else
        matlabbatch{4}.spm.spatial.coreg.estimate.other = {''};  
    end
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    % segmentation (struc)
    %matlabbatch{5}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{5}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Named File Selector: struc(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{5}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{5}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{5}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,1')};
    matlabbatch{5}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{5}.spm.spatial.preproc.tissue(1).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(1).warped = [1 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,2')};
    matlabbatch{5}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{5}.spm.spatial.preproc.tissue(2).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(2).warped = [1 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,3')};
    matlabbatch{5}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{5}.spm.spatial.preproc.tissue(3).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(3).warped = [1 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,4')};
    matlabbatch{5}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{5}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,5')};
    matlabbatch{5}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{5}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,6')};
    matlabbatch{5}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{5}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{5}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{5}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{5}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{5}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{5}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{5}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{5}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{5}.spm.spatial.preproc.warp.write = [0 1];
    % normalization (struct, lesion)
    matlabbatch{6}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = options.resliceStruct;
    matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w_';
    % normalization (func)
    matlabbatch{7}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    for i = 1:length(func)
        matlabbatch{7}.spm.spatial.normalise.write.subj.resample(i) = cfg_dep(strcat('Realign: Estimate & Reslice: Realigned Images (Sess ',num2str(i),')'), substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{i}, '.','cfiles'));
    end
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(i+1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = options.resliceFunc;
    matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{7}.spm.spatial.normalise.write.woptions.prefix = 'w';
    % smoothing (func)
    matlabbatch{8}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{8}.spm.spatial.smooth.fwhm = options.fwhm;
    matlabbatch{8}.spm.spatial.smooth.dtype = 0;
    matlabbatch{8}.spm.spatial.smooth.im = 0;
    matlabbatch{8}.spm.spatial.smooth.prefix = 's';

    % save preproc batch
    save(fullfile(outDir,'batch_preproc.mat'),'matlabbatch');
    % run preproc batch
    spm_jobman('run',matlabbatch);
    
    % move data to outDir
    % move normalized struc and individual TPMs
    outDirStruc = fullfile(outDir,'struc');
    [strucPth,strucName,strucExt] = fileparts(struc(1,:));
    myMoveFile(strucPth,outDirStruc,strcat('w_',strucName,strucExt));
    myMoveFile(strucPth,outDirStruc,strcat('wc1',strucName,strucExt));
    myMoveFile(strucPth,outDirStruc,strcat('wc2',strucName,strucExt));
    myMoveFile(strucPth,outDirStruc,strcat('wc3',strucName,strucExt));    
    % move normalized lesion
    if isfield(options,'lesion')
        [lesionPth,lesionName,lesionExt] = fileparts(options.lesion);
        movefile(fullfile(lesionPth,strcat('w_',lesionName,lesionExt)),fullfile(outDirStruc,'wc4_lesion.nii'));
        % delete masked structural image
        delete(struc);
    end
    % move meanEpi
    [func1Pth,~,~] = fileparts(func{1}{1});
    myMoveFile(func1Pth,fullfile(outDir,'func'),'wmean*.nii');
    % delete swmeanEpi to avoid later copying
    delete(fullfile(func1Pth,'swmean*.nii'));
    delete(fullfile(func1Pth,'mean*.nii'));
    for iCond = 1:length(cond)
        % move smoothed and unsmoothed func and rp_
        [funcPth,~,~] = fileparts(func{iCond}{1});
        myMoveFile(funcPth,fullfile(outDir,'func',cond{iCond}),'sw*.nii');
        myMoveFile(funcPth,fullfile(outDir,'func',cond{iCond}),'w*.nii');
        myMoveFile(funcPth,fullfile(outDir,'func',cond{iCond}),'rp_*.txt');
        if isfield(options,'sliceorder')
            delete(fullfile(funcPth,'a*.nii'));
        end
    end
    
    % plot movement and save plot
    for iCond = 1:length(cond)
        % move func
        rpDir = fullfile(outDir,'func',cond{iCond});
        d = dir(fullfile(rpDir,'rp_*.txt'));
        mov(iCond) = afxPlotMovement(fullfile(rpDir,d(1).name),fullfile(outDir,strcat('movement-',cond{iCond},'.png')));
    end
    
    % load previeus subject struct array
    subjectsFile = fullfile('data',sampleName,'subjects','subjects.mat');
    if exist(subjectsFile,'file')
        load(subjectsFile);
    else
        subjects = struct([]);
    end
    
    % prepare subject structure
    subjects(end+1).name  = subjectName;
    subjects(end).dir   = fullfile(sampleName,'subjects',subjectName);
    subjects(end).sampleName = sampleName;
    subjects(end).exclude = 0;
    subjects(end).excludeReason = {};
    subjects(end).masks = spm_select('FPList',outDirStruc,'^wc.*.nii');
    subjects(end).t1    = spm_select('FPList',outDirStruc,'^w_.*.nii');
    for iCond = 1:length(cond)
        outDirFunc = fullfile(outDir,'func',cond{iCond});
        subjects(end).conditions(iCond).name  = cond{iCond};
        subjects(end).conditions(iCond).movement = mov(iCond);
        subjects(end).conditions(iCond).rmsFD = mov(iCond).rmsFD;
        subjects(end).conditions(iCond).func  = spm_select('ExtFPList',outDirFunc,'^sw.*.nii',Inf);
        subjects(end).conditions(iCond).func2 = spm_select('ExtFPList',outDirFunc,'^w.*.nii',Inf);
        subjects(end).conditions(iCond).rp    = spm_select('FPList',outDirFunc,'^rp_.*.txt');
        % exclude if less than 5 minutes of rs fmri remain with FD < .5 mm
        if (length(mov(iCond).FD)-length(mov(iCond).fd5))*options.TR < 5*60
            subjects(end).exclude = 1;
            subjects(end).excludeReason{end+1} = ['Excessive motion in condition ' cond{iCond}];
        end
    end 
    subjects(end).preprocOptions = options;
    subjects(end).info = info;
    % save new subject struct array
    save(subjectsFile,'subjects');
end

function myMoveFile(source,dest,pattern)
	if ~exist(dest,'dir'), mkdir(dest); end
	movefile(fullfile(source,pattern),dest);
end