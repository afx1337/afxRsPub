function [conn,connAll,subsAll,conditions] = afxSecondlevelNetwork(firstlevelInfo,networkName,networkNodes)

fprintf('afxSecondlevelNetwork ...\n');
    % defaults
    if nargin < 1 || isempty(firstlevelInfo)
        firstlevelInfo = spm_select(1,'^firstlevel_info.mat$','Select firstlevel_info.mat',{},'results');
    end

    % load firstlevel information
    info = load(firstlevelInfo);
    info.subjects = info.subjects(~[info.subjects.exclude]); % get rid of excluded subjects
    
    % define default network
    if nargin < 2
        networkName = 'all';
        networkNodes = 1:length(info.rois);
    end
    nRoi = length(networkNodes);
    
    % get all conditions
    tmpCond = [info.subjects.conditions];
    conditions = unique({tmpCond.name});    
    
    % 2nd level directory and output files
    dirFirstlevel = fullfile('results',info.firstlevelDir,'firstlevel');
    dirSecondlevel = fullfile('results',info.firstlevelDir,'secondlevel','network',networkName);
    if ~exist(dirSecondlevel,'dir'), mkdir(dirSecondlevel); end
    fnameOutMat = fullfile(dirSecondlevel,'networkConnectivity.mat');
    fnameOutMatAll = fullfile(dirSecondlevel,'networkConnectivityAll.mat');
    fnameOutXls = fullfile(dirSecondlevel,'networkConnectivity.xls');
    
    % number of subjects and names
    nSubAll = length(info.subjects);
    subsAll = {info.subjects.name}';
     
    % number of rois and relevant indices of full correlation matrix
    indInter = logical(triu(repmat([ones(nRoi,1) zeros(nRoi,1)],1,nRoi/2),-1)-triu(repmat([ones(nRoi,1) zeros(nRoi,1)],1,nRoi/2),0));
    indIntraL = logical(repmat([repmat([1 0],1,nRoi/2); zeros(1,nRoi)],nRoi/2,1));
    indIntraR = logical(repmat([zeros(1,nRoi); repmat([0 1],1,nRoi/2)],nRoi/2,1));
    indHalf = logical(triu(ones(nRoi),0));
    indIntraL(indHalf) = 0;
    indIntraR(indHalf) = 0;
    indAll = indInter | indIntraL | indIntraR;
    
    for iCond = 1:length(conditions)
        fprintf('  - condition: %s\n',conditions{iCond});
        load(fullfile(dirFirstlevel,['cond_' conditions{iCond}],'network.mat'));
        
        z = z(networkNodes,networkNodes,:);
        roiLesioned = roiLesioned(:,networkNodes);
        roiNames = roiNames(networkNodes);
        subsCond = subjectNames;
        % calculate interhemispheric and left/right intrahemispheric connectivity
        % relies in the fact that left hemisphere rois have odd indices and are followed
        % dirctly by ther right hemisphere homologue

        % calculate mean intr-/interhemispheric connectivity per patient
        for iSubAll = 1:nSubAll
            iSubCond = find(strcmp(subsCond,info.subjects(iSubAll).name));
            if length(iSubCond) == 1
                % get subjects data
                tmp = z(:,:,iSubCond);
                % get partially lesioned rois
                lesioned2D = repmat(roiLesioned(iSubCond,:),nRoi,1);
                lesioned2D = max(lesioned2D,lesioned2D');
                indDist = lesioned2D == 0;
                indPeri = lesioned2D > 0 & lesioned2D < 1;
                % generate conn indices
                conn.dist.Inter(iSubAll,iCond)  = nanmean(tmp(indInter & indDist));
                conn.dist.IntraL(iSubAll,iCond) = nanmean(tmp(indIntraL & indDist));
                conn.dist.IntraR(iSubAll,iCond) = nanmean(tmp(indIntraR & indDist));
                conn.peri.Inter(iSubAll,iCond)  = nanmean(tmp(indInter & indPeri));
                conn.peri.IntraL(iSubAll,iCond) = nanmean(tmp(indIntraL & indPeri));
                conn.peri.IntraR(iSubAll,iCond) = nanmean(tmp(indIntraR & indPeri));
                connAll(iSubAll,:,iCond) = tmp(indAll);
            else
                conn.dist.Inter(iSubAll,iCond)  = NaN;
                conn.dist.IntraL(iSubAll,iCond) = NaN;
                conn.dist.IntraR(iSubAll,iCond) = NaN;
                conn.peri.Inter(iSubAll,iCond)  = NaN;
                conn.peri.IntraL(iSubAll,iCond) = NaN;
                conn.peri.IntraR(iSubAll,iCond) = NaN;
                connAll(iSubAll,:,iCond) = nan(1,nnz(indAll));
            end
        end
    end
    
    % save results
    [pth,~,~] = fileparts(fnameOutMat);
    if ~exist(pth,'dir'), mkdir(pth); end
    % mat
    subs = subsAll;
    save(fnameOutMat,'conn','subs','conditions','roiNames');
    save(fnameOutMatAll,'connAll','subs','conditions','indAll','roiNames','roiLesioned');
    % xls
    tmp = repmat({'inter_dist' 'inter_peri' 'intraL_dist' 'intraL_peri' 'intraR_dist' 'intraR_peri' },length(conditions),1);
    out = ['condition' repmat(conditions,1,6) ; ['name' tmp(:)'] ; subsAll num2cell([conn.dist.Inter conn.peri.Inter conn.dist.IntraL conn.peri.IntraL conn.dist.IntraR conn.peri.IntraR]) ];
    xlswrite(fnameOutXls,out);
    fprintf('afxSecondlevelNetwork ... done\n');
end