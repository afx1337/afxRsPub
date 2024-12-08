function afxALFF(y,brainMask,subjectMasks,denoisingOptions,dim,mat,outDir,subjectName)
    % calculation of fALFF
    % Zang et al. 2007 "Altered baseline brain activity in children with
    % ADHD revealed by resting-state functional MRI."
    % Zou et al. 2008 "An improved approach to detection of amplitude of
    % low-frequency fluctuation (ALFF) for resting-state fMRI: Fractional ALFF"
    fprintf('   Calculating whole brain ALFF ... ');
    
    % frequencies for ALFF
    f1 = [.01 .08];
    
    % detrending
    y = detrend(y);

    % from https://de.mathworks.com/help/matlab/ref/fft.html
    L = size(y,1);  % length of signal
    Fs = 1/denoisingOptions.TR; % fs = 1/TR
    % Compute the Fourier transform of the signal.
    yFft = fft(y,[],1)/L;
    % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(yFft);
    P1 = P2(1:floor(L/2+1),:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    % Define the frequency domain f and
    f = Fs*(0:(L/2))/L;
    % define the bins of interest for the amplitude analysis
    idx1 = f>=f1(1) & f<=f1(2);
    % init fALFF array
    ALFF = nan(1,size(y,2));
    % ALFF is defined as the average square root of the powers between .01
    % and .08 Hz
    % For fALFF this is divided by the same measure for the whole spectrum
    %fALFF(mask) = mean(sqrt(P1(idx1,mask))) ./ mean(sqrt(P1(idx2,mask)));
    ALFF(brainMask) = sum(P1(idx1,brainMask)); % adapted from rest toolbox
    % the resulting values are then transformed to z-scores (within the individual gm-mask)
    gmMask = subjectMasks(:,1) > .2 & ~isnan(ALFF)';
    ALFF(~isnan(ALFF)) = ALFF(~isnan(ALFF)) - mean(ALFF(gmMask));
    ALFF(~isnan(ALFF)) = ALFF(~isnan(ALFF))/std(ALFF(gmMask));
    
    % create output dir if necassary
    imgFname = fullfile(outDir,'ALFF',strcat(subjectName,'.nii'));
    [pth,~,~] = fileparts(imgFname);
    if ~exist(pth,'dir'), mkdir(pth); end
    % save fcunctional connectivity map
    afxVolumeWrite(imgFname,ALFF,dim,'int16',mat,'ALFF z-scores')
    fprintf('done\n')
end