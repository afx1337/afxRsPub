function afxPreprocHCP(options,cond,sampleName,subjectName,info,fix)
    
    if nargin < 5, info = struct([]); end
    if nargin < 6, fix = false; end

    outDir = fullfile(pwd,'data',sampleName,'subjects',subjectName);
    if fix
        movefile(outDir,strcat(outDir,'_broken'))
    else
        if exist(fullfile(outDir,'func'),'dir')
            disp(['Skipping subject >' subjectName '< (preprocessing has already been done)']);
            return;
        end
    end
    if ~exist(outDir,'dir'), mkdir(outDir); end
    
     % create temporary directory
    dTemp = tempname(); mkdir(dTemp);
    
    % prepare subject structure
    thisSubject.name  = subjectName;
    thisSubject.dir   = fullfile(sampleName,'subjects',subjectName);
    thisSubject.sampleName = sampleName;
    thisSubject.exclude = 0;
    thisSubject.excludeReason = {};
    thisSubject.masks = spm_select('FPList','masks','^wc.*.nii');
    thisSubject.t1 = '';
    thisSubject.preprocOptions = options;
    thisSubject.info = info;
    
    for iCond = 1:length(cond) % condition 1 is lr, 2 is rl
        outDirFunc = fullfile(outDir,'func',cond(iCond).name);
        mkdir(outDirFunc);
        
        % unzip if nii is gziped
        if strcmp(cond(iCond).func(end-6:end),'.nii.gz')
            gunzip(cond(iCond).func,dTemp);
            [~,f,~] = fileparts(cond(iCond).func);
            cond(iCond).func = fullfile(dTemp,f);
        end
        % load rsfMRI images
        Nfunc = nifti(cond(iCond).func);
        % get image properties
        nVox = prod(Nfunc.dat.dim(1:3));
        Y = nan(Nfunc.dat.dim(4),nVox);
        parfor j = 1:Nfunc.dat.dim(4)
        	Y(j,:) = reshape(Nfunc.dat(:,:,:,j),1,[]);
        end % ~50 seconds per rs session
        dim = Nfunc.dat.dim(1:3);
        mat = Nfunc.mat;
        [R,C,P] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
        RCP     = [R(:)';C(:)';P(:)';ones(1,numel(R))];
        XYZmm   = mat(1:3,:)*RCP;
        XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
        clear Nfunc R C P RCP dat
        % clear temporary directory
        delete(fullfile(dTemp,'*'));
        
        % convolve functional images with gaussian smoothing kernel
        Ys = afxSmoooth(Y,options.fwhm,dim,mat); % ~85s
        
        % save mean rest
        afxVolumeWrite(fullfile(outDirFunc,'meanwrest.nii'),mean(Y),dim(1:3),'int16',mat,'HCP preprocessed rest',true);
        
        % create 1d-niftis using brainmask.nii to save disk-space and to
        % speed up loading processes
        [brainMask,~] = afxLoadMasks([],fullfile('masks','devMaskHCP.nii'),XYZmm);
        Y = Y(:,brainMask);
        Ys = Ys(:,brainMask);
        if max(Y(:)) > intmax('uint16') || max(Y(:)) < intmin('uint16')
            error('Data out of UINT16 range');
        end
        Y = uint16(Y); % rounding errors are relatively small
        Ys = uint16(Ys);
       
        % prepare motion parametes to rp_-style
        rp = dlmread(cond(iCond).mov);
        rp = rp(:,1:6); % trns_x trans_y trnas_z rot_x rot_y rot_z (mm/deg)
        rp(:,4:6) = rp(:,4:6)*pi()/180; % rot deg->rad
        % save rp_*.txt
        rpFile = fullfile(outDirFunc,'rp_rest.txt');
        dlmwrite(rpFile,rp,'\t');
        
        % plot movement and save plot
        mov(iCond) = afxPlotMovement(rpFile,fullfile(outDir,strcat('movement-',cond(iCond).name,'.png')));

        % save files to data dir
        wrestFile = fullfile(outDirFunc,'wrest.mat');
        swrestFile = fullfile(outDirFunc,'swrest.mat');
        save(wrestFile,'Y','mat','dim','brainMask','-v7.3','-nocompression');
        Y = Ys;
        save(swrestFile,'Y','mat','dim','brainMask','-v7.3','-nocompression');
        
        thisSubject.conditions(iCond).name  = cond(iCond).name;
        thisSubject.conditions(iCond).movement = mov(iCond);
        thisSubject.conditions(iCond).rmsFD = mov(iCond).rmsFD;
        thisSubject.conditions(iCond).func  = swrestFile;
        thisSubject.conditions(iCond).func2 = wrestFile;
        thisSubject.conditions(iCond).rp    = rpFile;
        % exclude if less than 10 minutes of rs fmri remain with FD < .5 mm
        if (length(mov(iCond).FD)-length(mov(iCond).fd5))*options.TR < 10*60
            thisSubject.exclude = 1;
            thisSubject.excludeReason{end+1} = ['Excessive motion in condition ' cond(iCond).name];
        end
    end

    % load previeus subject struct array and append
    if ~fix
        subjectsFile = fullfile('data',sampleName,'subjects','subjects.mat');
        if exist(subjectsFile,'file')
            load(subjectsFile);
            subjects(end+1) = thisSubject;
        else
            % ... or create
            subjects = thisSubject;
        end

        % save subject struct array
        save(subjectsFile,'subjects');
    end
    
    % clear temporary directory
    rmdir(dTemp);
end