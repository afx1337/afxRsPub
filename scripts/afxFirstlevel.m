function projectDir = afxFirstlevel(subjects,denoisingOptions,rois,projectName,analyses)

    % default analyses
    if nargin < 5
        analyses = {'fc_wholebrain','fc_network'};
    end    

    % load subjects, options and rois if necassary
    if ischar(subjects), load(subjects); end
    if ischar(denoisingOptions), load(denoisingOptions); end
    if ischar(rois), load(rois); end
    [rois.ind] = deal([]);
    
    % update paths in case data has moved
    subjects = afxUpdatePaths(subjects);
    
    % write firstlevel to results\sampleName\projectName
    % (sampleName is taken from the first subject)
    firstlevelDir = fullfile(subjects(1).sampleName,projectName);
    projectDir = fullfile('results',firstlevelDir);

    % calculate firstlevel for all subjects and conditions
    includedSubjects = find(~[subjects.exclude]);
    for iSubject = includedSubjects
        for iCond = 1:length(subjects(iSubject).conditions)
            denoisingOptions.TR = subjects(iSubject).preprocOptions.TR;
            fprintf('Subject >%s< (%i/%i) condition >%s< ...\n',subjects(iSubject).name,iSubject,length(subjects),subjects(iSubject).conditions(iCond).name);
            % add empty func2-fieled to ensure backwards compatibility
            if ~isfield(subjects(iSubject).conditions(iCond),'func2'), subjects(iSubject).conditions(iCond).func2 = []; end
            rois = afxConn(...
                struct('func',subjects(iSubject).conditions(iCond).func,'func2',subjects(iSubject).conditions(iCond).func2),...
                subjects(iSubject).masks,...
                subjects(iSubject).conditions(iCond).rp,...
                denoisingOptions,...
                rois,...
                firstlevelDir,...
                subjects(iSubject).conditions(iCond).name,...
                subjects(iSubject).name,...
                analyses...
            );
            fprintf('\n');
        end
    end
    
    % save project information
    meta = afxMetaData();
    rois = rmfield(rois,'ind');
    save(fullfile(projectDir,'firstlevel_info.mat'),'subjects','rois','denoisingOptions','projectName','firstlevelDir','analyses','meta');
    
    if any(contains(analyses,'piniroi'))
        afxSecondlevelPiniROIs(fullfile(projectDir,'firstlevel_info.mat'))
    end
end
