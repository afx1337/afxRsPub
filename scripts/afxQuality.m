function afxQuality(y,subjectMasks,dim,mat,outDir,subject)
    % calculation of quality checks
    fprintf('   Calculating quality checks ... ')
    
    % calculate percent signal change (psc)
    yMean = mean(y);
    y = y./yMean*100-100;
    % calculate root mean square of psc
    yRms = sqrt(mean(y.^2));
    % get mean rms psc within grey matter
    meanRms = nanmean(yRms(subjectMasks(:,1) > .5));
    
     % create output dir if necassary
     imgFname = fullfile(outDir,'quality',strcat(subject,'.nii'));
     [pth,~,~] = fileparts(imgFname);
     if ~exist(pth,'dir'), mkdir(pth); end
    
     % save rms psc map
     afxVolumeWrite(imgFname,yRms,dim,'int16',mat,'rms of psc');
     
%      % generate image
%      pngFname = fullfile(outDir,'quality',strcat(subject,'.png'));
%      f = afxPlotSections(imgFname,-55:10:80,'',[1 Inf],[1 5],1);
%      f.PaperUnits = 'inches';
%      f.PaperPosition = [0 0 1367 1024]./100;
%      figure(f);
%      print(pngFname,'-dpng','-r100');
%      close(f);

    % filenames for output
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fnameOutTxt = fullfile(outDir,'quality.txt');
    fnameOutMat = fullfile(outDir,'quality.mat');

    % try to load previus data
    if exist(fnameOutMat,'file')
        load(fnameOutMat);
        % check if subject is already present
        ind = strcmp(subjectNames,subject);
        if any(ind)
            % replace
            meanGmRmsPsc(ind) = meanRms;
        else
            % append
            subjectNames{end+1} = subject;
            meanGmRmsPsc(end+1) = meanRms;
        end
    else
        % new
        subjectNames = {subject};
        meanGmRmsPsc(1) = meanRms;
    end

   afxCsvWrite(fnameOutTxt,[['subject'; subjectNames'] ['meanGmRmsPsc'; num2cell(meanGmRmsPsc')]]);
    save(fnameOutMat,'meanGmRmsPsc','subjectNames');
    
    fprintf('done\n')
end