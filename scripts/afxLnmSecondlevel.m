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
    
    % load underlay for rendering
    [templateVol,~,templateDim,templateMat] = afxLoadFunc(char(fullfile('templates','MNI152_T1_0.5mm_masked.nii')));
    templateVol = reshape(templateVol, templateDim);
    opts.intensityWindow = [0 255];

    % rendering defaults
    sliceSpec.orientation = 'axial';
    sliceList = [40 20 0 -30];
    overlays(1).thr = [tCrit Inf];
    overlays(1).mode = 'gradient';
    overlays(1).colormap = 'hot';

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
            overlays(1).fname = fT;
            afxRenderSections(templateVol, templateMat, sliceSpec, sliceList, overlays, opts, fTPng);
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
end
