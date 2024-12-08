function ztabledat = afxConnNetwork(yRoi,lesioned,rois,outPrefix,subject)
    % calculation of functional connectivity
    fprintf('   Calculating network functional connectivity ... ')
    
    % calculate correlation coefficients and apply fisher transformation
    rhoFisher = atanh(corrcoef(yRoi));
    % write data to ztable if <= 10 rois
    if size(rhoFisher,1) <= 10
        n = length(rhoFisher);
        ztabledat = {subject}; header = {'subject'};
        for iCol = 1:(n-1)
            for iRow = (iCol+1):n
                header{end+1} = [ rois(iCol).name '<->' rois(iRow).name ];
                ztabledat{end+1} = rhoFisher(iRow,iCol);
            end
        end
    else
        ztable = 'No z table for more than 10 ROIs';
    end
    
    % filenames for output
    fnameOutTxt = fullfile(outPrefix,'network.txt');
    fnameOutMat = fullfile(outPrefix,'network.mat');

    % try to load previus data ztable
    if exist(fnameOutMat,'file')
        load(fnameOutMat);
        % check if subject is already present
        ind = strcmp(subjectNames,subject);
        if any(ind)
            % replace
            z(:,:,ind) = rhoFisher;
            roiLesioned(ind,:) = lesioned;
            if size(rhoFisher,1) <= 10, ztable(1+find(ind),:) = ztabledat; end
        else
            % append
            subjectNames{end+1} = subject;
            z(:,:,end+1) = rhoFisher;
            roiLesioned(end+1,:) = lesioned;
            if size(rhoFisher,1) <= 10, ztable = [ztable; ztabledat]; end
        end
    else
        % new
        subjectNames = {subject};
        z(:,:,1) = rhoFisher;
        roiLesioned(1,:) = lesioned;
        if size(rhoFisher,1) <= 10, ztable = [header; ztabledat]; end
    end
    
    roiNames = {rois.name};
    
    % write to txt and mat file
    [pth,~,~] = fileparts(fnameOutMat);
    if ~exist(pth,'dir'), mkdir(pth); end
    
    if size(rhoFisher,1) <= 10, afxCsvWrite(fnameOutTxt,ztable); end
    save(fnameOutMat,'ztable','z','roiLesioned','roiNames','subjectNames');
    
    fprintf('done\n')
end