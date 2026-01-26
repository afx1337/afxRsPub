function [rois, lesioned, yRoi, yRoiCC] = afxRoiTimeseries(rois,y,subjectMasks,brainMask,XYZmm,denoisingOptions)

    fprintf('   ROI timeseries extraction ...')
    yRoi = nan(size(y,1),length(rois));
    lesioned = zeros(length(rois),1);

    % initialize caches (avoid repeated resampling of atlas or gm file for
    % multiple ROIs
    atlasCache = struct('file',[],'dat',[]);
    gmCache = struct('file',[],'dat',[]);
    
    % iterate over rois
    for iRoi = 1:length(rois)
        % get ROI voxel indices
        if isempty(rois(iRoi).ind)
            switch rois(iRoi).type
                case 'sphere'
                    xyz = rois(iRoi).coords;
                    D = sqrt((XYZmm(1,:)-xyz(1)).^2+(XYZmm(2,:)-xyz(2)).^2+(XYZmm(3,:)-xyz(3)).^2);
                    rois(iRoi).ind = (D < rois(iRoi).radius)';
                    clear D;
                case 'image'
                    dat = afxVolumeResample(rois(iRoi).file,XYZmm,0);
                    rois(iRoi).ind = dat > 0.5;
                    clear dat;
                case 'atlas'
                    atlasCache = getResampledData(atlasCache, rois(iRoi).file, XYZmm);
                    rois(iRoi).ind = atlasCache.dat > rois(iRoi).pick-.5 & atlasCache.dat < rois(iRoi).pick+.5 ;
                otherwise
                    error('unknown roi type: %s',rois(iRoi).type);
            end
        end

        % check if roi is lesioned
        if size(subjectMasks,2) > 3
            lesioned(iRoi) = nnz(rois(iRoi).ind & subjectMasks(:,4)>=.5) / nnz(rois(iRoi).ind);
        end

        % perform brain (and implicitly lesion) masking (just to be sure); proceed with copy of ind
        ind = rois(iRoi).ind & brainMask;
        
        
        % save cross-correlation matrix before gm masking (Pini et al. 2022)
        if nargout > 3
            yRoiCC(iRoi).name = rois(iRoi).name;
            cc = atanh(corr(y(:,ind)));
            cc(eye(size(cc),'logical')) = NaN;
            yRoiCC(iRoi).cc = nanmean(cc,2);
            yRoiCC(iRoi).ind = ind;
            continue;
        end
        
        % perform gray matter masking
        switch rois(iRoi).gmMask
            case 0 % gm no masking
            case 'split' % median split
                ind = ind & subjectMasks(:,1) >= median(subjectMasks(ind,1));
            otherwise
                if ischar(rois(iRoi).gmMask)
                    if ~exist(rois(iRoi).gmMask, 'file')
                        error('gmMask file not found: %s',rois(iRoi).gmMask);
                    end
                    gmCache = getResampledData(gmCache, rois(iRoi).gmMask, XYZmm);
                    ind = ind & (gmCache.dat ~= 0);
                elseif isnumeric(rois(iRoi).gmMask) && rois(iRoi).gmMask > 0 && rois(iRoi).gmMask < 1
                    ind = ind & subjectMasks(:,1) > rois(iRoi).gmMask;
                else
                    error('unknown gmMask option: %s', mat2str(rois(iRoi).gmMask));
                end
        end
        
        % extract timeseries
        ts = y(:,ind);
        
        if any(isnan(ts(:)))
            ts = zeros(size(y,1),1);
            warning(['ROI ' rois(iRoi).name ' contains NaNs. Proceeding with zero timeseries ...']);
        elseif isempty(ts) || sum(var(ts)) == 0
            ts = zeros(size(y,1),1);
            warning(['ROI ' rois(iRoi).name ' is empty or has no variance. Proceeding with zero timeseries ...']);
        else
            if denoisingOptions.unsmoothed
                % mean over voxels for unsmoothed data
                ts = mean(ts,2);
            else
                % first eigenvariate for smooth data
                ts = afxEigenvariate(ts);
            end
        end
        % save
        yRoi(:,iRoi) = ts;
    end
    
    fprintf(' done\n')
end

function cache = getResampledData(cache, file, XYZmm)
    % Funktion zur Vermeidung mehrfacher Resampling-Operationen auf das gleiche File
    if ~strcmp(cache.file, file)
        cache.dat = afxVolumeResample(file, XYZmm, 0); % Resample
        cache.file = file; % Dateinamen speichern
    end
end
