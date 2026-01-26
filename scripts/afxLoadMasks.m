function [brainMask,subjectMasks] = afxLoadMasks(fmasks,fbrainMask,XYZmm)

    fprintf('   Loading masks ... ')
    % load (individual) tpms and resample to func space (trilinear interpolation)
    % if necessary
    subjectMasks = zeros(size(XYZmm,2),size(fmasks,1),'single');
    for  i = 1:size(fmasks,1)
        subjectMasks(:,i) = afxVolumeResample(fmasks(i,:),XYZmm,1);
    end
    % load brainmask and resample to func space (nearest neighbur) if
    % necessary and binarize
    if ~isempty(fbrainMask)
        brainMask = afxVolumeResample(fbrainMask,XYZmm,0)>.5;
    else
        brainMask = [];
    end
    % subtract lesion from brain mask
    if size(subjectMasks,2) > 3
        brainMask = brainMask & subjectMasks(:,4)'<.5;
    end
    fprintf('done\n');
end