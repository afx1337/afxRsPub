function afxConnHemodynamicLag(y,brainMask,subjectMasks,rpFile,options,dim,mat,outDir,subjectName)
    % calculation of hemodynamic lag
    % See Lv et al. (2012) and  Siegel et al. (2015)
    fprintf('   Calculating hemodynamic lag ...\n')

    % calculate global signal (from GM, see Siegel et al., 2015)
    indGSM = subjectMasks(:,1) > .95; %| subjectMasks(:,2)>options.threshWM;
    gs(:,1) = mean(y(:,indGSM).*repmat(subjectMasks(indGSM,2),1,size(y,1))',2);
    
    % mean centering of global signal
    gs(:,1) = gs(:,1) - mean(gs(:,1));
    
    delays = -4:4;
    
    % get temporal mask (FD <= .5 mm)
    rp = dlmread(rpFile);
    rpp = [zeros(1,size(rp,2)); diff(rp)];
    FD = sum(abs([rpp(:,1:3) rpp(:,4:6)*50]),2);
    temporalMask = FD <= .5;
    for iDelay = 1:length(delays)
        if delays(iDelay) < 0
            a = -delays(iDelay);
            b = 0;
        else
            a = 0;
            b = delays(iDelay);
        end
        tm{iDelay} = temporalMask(a+1:end-b) & temporalMask(b+1:end-a);
    end
    
    tr = options.TR;
    lag = nan(1,size(y,2));
    parfor i = 1:size(y,2)
        cc = nan(1,length(delays));
        if brainMask(i)
            for iDelay = 1:length(delays)
                if delays(iDelay) < 0
                    a = -delays(iDelay);
                    b = 0;
                else
                    a = 0;
                    b = delays(iDelay);
                end
                % temporal shift
                s = y(a+1:end-b,i);
                g = gs(b+1:end-a);
                s = s(tm{iDelay});
                g = g(tm{iDelay});
                % lagged cross-correlation, see Siegel et al., 2015
                % mean-centering???
                cc(iDelay) = 1/length(g) * g'*s./ (std(g) * std(s));
            end
            % get maximum cross-correlation
            [~,ind] = max(cc);
            if ind > 1 && ind < length(delays)
                % fit parabolic function
                ccmax = afxParabolicIntMax(delays(ind-1:ind+1),cc(ind-1:ind+1))*tr;
            else
                ccmax = NaN; %delays(ind)*tr;
            end
            lag(i) = -ccmax;
        end
    end
    
    % reshape and mirror lesion
    lesion = reshape(subjectMasks(:,4)>.5,dim);
    lesion = lesion | flip(lesion,1);
   
    lag3d = reshape(lag,dim);
    % ignore lag within lesion and mirrored lesion
    lag3d(lesion) = NaN;
    %lag3d(~indGSM) = NaN;
    % smoothing
    %lag3d = smooth3(lag3d,'gaussian',1);
    %lag = lag3d(:);

    % calculate laterality
    lat = mat(1,1)/abs(mat(1,1));
    mid1 = floor(dim(1)/2);
    mid2 = ceil(dim(1)/2);
    lh = lag3d(1:mid1,:,:);
    rh = lag3d(mid2:end,:,:);
    clear lag3d;
    lh = mean(lh(~isnan(lh)));
    rh = mean(rh(~isnan(rh)));

    LL = lat*(lh - rh);  % = average lag difference in seconds
    %fprintf('LI: %f; LL: %f',LI,LL);

    % create output dir if necassary
    imgFname = fullfile(outDir,'hemo_lag',strcat(subjectName,'.nii'));
    [pth,~,~] = fileparts(imgFname);
    if ~exist(pth,'dir'), mkdir(pth); end
    % save lag map
    afxVolumeWrite(imgFname,lag,dim,'int16',mat,'hemodynamic lag in s');

    % write to file
    fnameOutTxt = fullfile(outDir,'hemo_lag.txt');
    fnameOutMat = fullfile(outDir,'hemo_lag.mat');
    header = {'subject', 'lag laterality (s)'};
    % try to load previus data ztable
    if exist(fnameOutMat,'file')
        load(fnameOutMat);
        allSubs = lag_laterality(2:end,1);
        % check if subject is already in the ztable
        ind = strcmp(allSubs,subjectName);
        if ~isempty(find(ind,1))
            % replace
            lag_laterality(1+find(ind),:) = { subjectName LL };
        else
            % appand
            lag_laterality = [lag_laterality; { subjectName LL } ];
        end
    else
        % new table
        lag_laterality = [header; { subjectName LL } ];
    end
    
    % write to txt and mat file
    [pth,~,~] = fileparts(fnameOutMat);
    if ~exist(pth,'dir'), mkdir(pth); end
    afxCsvWrite(fnameOutTxt,lag_laterality);
    save(fnameOutMat,'lag_laterality');
    
    fprintf('done\n')
end

function x_extr = afxParabolicIntMax(x,y)
    X = [(x.^2)' x' ones(3,1)]; %Matrix of x-terms
    constants = X\(y)';
    x_extr = -.5*(constants(2)/constants(1)); %From the first derivative
end
