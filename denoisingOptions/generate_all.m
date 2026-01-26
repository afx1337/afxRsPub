clear
denoisingOptions = struct;
denoisingOptions.regressRp = false;
denoisingOptions.regressTs = struct([]);
denoisingOptions.filter = [0 Inf];
denoisingOptions.threshFD = Inf;
denoisingOptions.unsmoothed = false;
save('unsmoothed_quality.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .95;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .95;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = false;
save('smoothed_95_WMCSF.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .95;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .95;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_95_WMCSF.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .95;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .95;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.regressTs(3).tpm = [1 1 0]; % GS
denoisingOptions.regressTs(3).thresh = .95;
denoisingOptions.regressTs(3).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = false;
save('smoothed_95_GSR.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .95;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .95;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.regressTs(3).tpm = [1 1 0]; % GS
denoisingOptions.regressTs(3).thresh = .95;
denoisingOptions.regressTs(3).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_95_GSR.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM PCA denoising
denoisingOptions.regressTs(1).thresh = .99;
denoisingOptions.regressTs(1).pca = 5;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF PCA denoising
denoisingOptions.regressTs(2).thresh = .99;
denoisingOptions.regressTs(2).pca = 5;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = Inf; % no motion scrubbing for CompCor
denoisingOptions.unsmoothed = false;
save('smoothed_99_PCA.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM PCA denoising
denoisingOptions.regressTs(1).thresh = .99;
denoisingOptions.regressTs(1).pca = 5;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF PCA denoising
denoisingOptions.regressTs(2).thresh = .99;
denoisingOptions.regressTs(2).pca = 5;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = Inf; % no motion scrubbing for CompCor
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_99_PCA.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = false;
denoisingOptions.regressTs = struct([]);
denoisingOptions.filter = [0 Inf];
denoisingOptions.threshFD = Inf;
denoisingOptions.unsmoothed = false;
save('smoothed_fALFF.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .95;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .95;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.filter = [0 Inf];
denoisingOptions.threshFD = Inf;
denoisingOptions.unsmoothed = false;
save('smoothed_fALFF_95_WMCSF.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .80;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .85;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = false;
save('smoothed_HCP_WMCSF.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .80;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .85;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_HCP_WMCSF.mat');


clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .80;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .85;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.regressTs(3).tpm = [1 1 0]; % GS
denoisingOptions.regressTs(3).thresh = .5;
denoisingOptions.regressTs(3).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = false;
save('smoothed_HCP_GSR.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM
denoisingOptions.regressTs(1).thresh = .80;
denoisingOptions.regressTs(1).pca = 0;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF
denoisingOptions.regressTs(2).thresh = .85;
denoisingOptions.regressTs(2).pca = 0;
denoisingOptions.regressTs(3).tpm = [1 1 0]; % GS
denoisingOptions.regressTs(3).thresh = .5;
denoisingOptions.regressTs(3).pca = 0;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = 0.5;
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_HCP_GSR.mat');

clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM PCA denoising
denoisingOptions.regressTs(1).thresh = .90;
denoisingOptions.regressTs(1).pca = 5;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF PCA denoising
denoisingOptions.regressTs(2).thresh = .90;
denoisingOptions.regressTs(2).pca = 5;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = Inf; % no motion scrubbing for CompCor
denoisingOptions.unsmoothed = false;
save('smoothed_HCP_PCA.mat');


clear
denoisingOptions = struct;
denoisingOptions.regressRp = true;
denoisingOptions.regressTs(1).tpm = [0 1 0]; % WM PCA denoising
denoisingOptions.regressTs(1).thresh = .9;
denoisingOptions.regressTs(1).pca = 5;
denoisingOptions.regressTs(2).tpm = [0 0 1]; % CSF PCA denoising
denoisingOptions.regressTs(2).thresh = .9;
denoisingOptions.regressTs(2).pca = 5;
denoisingOptions.filter = [0.01 0.08];
denoisingOptions.threshFD = Inf; % no motion scrubbing for CompCor
denoisingOptions.unsmoothed = true;
denoisingOptions.sFWHM = 5; % mm
save('unsmoothed_HCP_PCA.mat');
