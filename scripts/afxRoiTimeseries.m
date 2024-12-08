function [yRoi, lesioned, yRoiCC] = afxRoiTimeseries(rois,y,subjectMasks,brainMask,XYZmm,denoisingOptions)

    fprintf('   ROI timeseries extraction ...')
    yRoi = nan(size(y,1),length(rois));
    lesioned = zeros(length(rois),1);

    % initialize atlas cache (avoid repeated resampling of atlas file for
    % multiple ROIs
    atlasCache = struct('file',[],'dat',[]);
    
    % iterate over rois
    for iRoi = 1:length(rois)
        % get ROI voxel indices
        curRoi = rois(iRoi);
        switch curRoi.type
            case 'sphere'
                xyz = curRoi.coords;
                D = sqrt((XYZmm(1,:)-xyz(1)).^2+(XYZmm(2,:)-xyz(2)).^2+(XYZmm(3,:)-xyz(3)).^2);
                ind = (D < curRoi.radius)';
                clear D;
            case 'image'
                dat = afxVolumeResample(curRoi.file,XYZmm,0);
                ind = dat > 0.5;
                clear dat;
            case 'atlas'
                if ~strcmp(atlasCache.file,curRoi.file)
                    atlasCache.dat = afxVolumeResample(curRoi.file,XYZmm,0);
                    atlasCache.file = curRoi.file;
                end
                ind = atlasCache.dat > curRoi.pick-.5 & atlasCache.dat < curRoi.pick+.5 ;
                clear dat;
            otherwise
                error('unknown roi type: %s',curRoi.type);
        end

        % check if roi is lesioned
        if size(subjectMasks,2) > 3
            lesioned(iRoi) = nnz(ind & subjectMasks(:,4)>=.5) / nnz(ind);
        end

        % perform brain (and implicitly lesion) masking (just to be sure)
        ind = ind & brainMask;
        
        % save cross-correlation matrix before gm masking (Pini et al. 2022)
        if nargout > 2
            yRoiCC(iRoi).name = curRoi.name;
            cc = atanh(corr(y(:,ind)));
            cc(eye(size(cc),'logical')) = NaN;
            yRoiCC(iRoi).cc = nanmean(cc,2);
            yRoiCC(iRoi).ind = ind;
        end
        
        % perform gray matter masking
        switch curRoi.gmMask
            case 0 % gm no masking
                ind = ind;
            case 'split' % median split
                ind = ind & subjectMasks(:,1)>=median(subjectMasks(ind,1));
            otherwise
                if isnumeric(curRoi.gmMask) && curRoi.gmMask > 0 && curRoi.gmMask < 1
                    ind = ind & subjectMasks(:,1)>curRoi.gmMask;
                else
                    error('unknown gmMask option: %s',curRoi.gmMask);
                end
        end
        
        % extract timeseries
        ts = y(:,ind);
        
        if isempty(ts)
            % flat timeseries if no voxels were selected, print warning
            ts = zeros(size(y,1),1);
            warning(['ROI ' curRoi.name ' is empty. (Maybe due to lesion or gm masking.) Proceeding with zero timeseries ...']);
        elseif sum(var(ts)) == 0
            % flat timeseries if variance is zero in all voxels. print warning
            ts = zeros(size(y,1),1);
            warning(['ROI ' curRoi.name ' contains no variance. (Maybe outside the brain.) Proceeding with zero timeseries ...']);
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