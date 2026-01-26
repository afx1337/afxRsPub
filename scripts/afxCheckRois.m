function [rois, delInd] = afxCheckRois(rois,subjects)

    fprintf('   Check for empty ROIs ...')

     gmCache = struct('file',[],'dat',[]);
    
    % get functional space
    [~,XYZmm,~,~] = afxLoadFunc(subjects(1).conditions(1).func(1,:));
    % get brain mask
    [brainMask,~] = afxLoadMasks([],fullfile('masks','brainmask.nii'),XYZmm);
    
    % load (all) subjects gm mask
    for iSubject = find(~[subjects.exclude])
        [~,subjectMasks{iSubject}] = afxLoadMasks(subjects(iSubject).masks(1,:),[],XYZmm);
    end
    
    delInd = false(1,length(rois));
    
    % check all rois (reversed order to allow direct deletion in struct array)
    for iRoi = length(rois):-1:1
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
            otherwise
                error('unknown roi type: %s',curRoi.type);
        end
        
        % perform brain masking (just to be sure)
        ind = ind & brainMask;
        
        % check if ROI is empty before gm masking
        if sum(ind) == 0
            fprintf('      > delete empty ROI %s\n',rois(iRoi).name);
            rois(iRoi) = [];
            delInd(iRoi) = true;
        elseif curRoi.gmMask ~= 0
            % perform gray matter masking for each mask
            for iSubject = find(~[subjects.exclude])
                switch curRoi.gmMask
                    case 0 % gm no masking
                        ind = ind;
                    case 'split' % median split
                        ind = ind & subjectMasks{iSubject}(:,1)>median(subjectMasks{iSubject}(ind,1));
                    otherwise
                        if ischar(curRoi.gmMask)
                            if ~exist(curRoi.gmMask, 'file')
                                error('gmMask file not found: %s',curRoi.gmMask);
                            end
                            gmCache = getResampledData(gmCache, curRoi.gmMask, XYZmm);
                            ind = ind & (gmCache.dat ~= 0);
                        elseif isnumeric(curRoi.gmMask) && curRoi.gmMask > 0 && curRoi.gmMask < 1
                            ind = ind & subjectMasks{iSubject}(:,1) > curRoi.gmMask;
                        else
                            error('unknown gmMask option: %s', mat2str(curRoi.gmMask));
                        end

                end
                % check if ROI is empty after masking
                if sum(ind) == 0
                    fprintf('      > delete empty ROI %s\n',rois(iRoi).name);
                    rois(iRoi) = [];
                    delInd(iRoi) = true;
                    break;
                end
            end
        end
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
