function subjects = afxUpdatePaths(subjects)
    baseDir = fullfile(pwd(),'data');
    
    % return, if t1 of first subject exists
    %dirCond = fullfile('func',subjects(1).conditions(1).name);
    %tmp = afxUpdate(subjects(1).conditions(1).func(1,:),dirCond,fullfile(baseDir,subjects(1).dir));
    if exist(subjects(1).conditions(1).func,'file') || exist(subjects(1).t1,'file'), return; end
    
    fprintf('Preprocessed data has moved. Updating all paths ... ');
    for iSub = 1:length(subjects)
        baseDirSub = fullfile(baseDir,afxFixPath(subjects(iSub).dir));
        % update t1
        subjects(iSub).t1 = afxUpdate(subjects(iSub).t1,'struc',baseDirSub);
        % update masks
        clear tmp;
        [p,~,~] = fileparts(afxFixPath(subjects(iSub).masks(1,:)));
        [~,p,~] = fileparts(p);
        if strcmp(p,'masks')
            for iMask = 1:size(subjects(iSub).masks,1)
                tmp{iMask} = afxUpdate(subjects(iSub).masks(iMask,:),'masks',pwd);
            end
        else
            for iMask = 1:size(subjects(iSub).masks,1)
                tmp{iMask} = afxUpdate(subjects(iSub).masks(iMask,:),'struc',baseDirSub);
            end
        end
        subjects(iSub).masks = cell2mat(tmp');
        % update conditions
        for iCond = 1:length(subjects(iSub).conditions)
            dirCond = fullfile('func',subjects(iSub).conditions(iCond).name);
            % update rp
            subjects(iSub).conditions(iCond).rp = afxUpdate(subjects(iSub).conditions(iCond).rp,dirCond,baseDirSub);
            % update func
            clear tmp;
            for iFunc = 1:size(subjects(iSub).conditions(iCond).func,1)
                  tmp{iFunc} = afxUpdate(subjects(iSub).conditions(iCond).func(iFunc,:),dirCond,baseDirSub);
            end
            subjects(iSub).conditions(iCond).func = cell2mat(tmp');
            % update func2
            clear tmp;
            for iFunc = 1:size(subjects(iSub).conditions(iCond).func2,1)
                 tmp{iFunc} = afxUpdate(subjects(iSub).conditions(iCond).func2(iFunc,:),dirCond,baseDirSub);
            end
            subjects(iSub).conditions(iCond).func2 = cell2mat(tmp');
        end
    end
    fprintf('done.\n');
end

function fileNew = afxUpdate(fileOld,subfolder,baseDir)
    fileOld = afxFixPath(fileOld);
    fileOld = afxFixPath(fileOld);
	[~,f,ext] = fileparts(fileOld);
    fileNew = fullfile(baseDir,subfolder,[f ext]);
end

function pth = afxFixPath(pth)
	pth = strrep(pth,'\',filesep);
    pth = strrep(pth,'/',filesep);
end