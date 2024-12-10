function afxConnVMHC(y,brainMask,dim,mat,outDir,subjectName)
    % calculation of VMHC
    fprintf('   Calculating VMHC ...\n')
    z = nan(1,size(y,2));
    parfor i = 1:size(y,2)
        % a = a - mean(a); b = b - mean(b)
        % rho = a*b / ( norm(a)*norm(b) )
        % note, that mean centering is not necessary for voxel timeseries
        % due to highpass filtering
        if brainMask(i)
            % get homotopic voxel index
            [ix1,iy,iz] = ind2sub(dim,i);
            ix2 = dim(1)-ix1+1;
            i2 = sub2ind(dim,ix2,iy,iz);
            % calculate conn
            if brainMask(i2)
                z(i) =  y(:,i2)' * y(:,i) ./ ( norm(y(:,i)) * norm(y(:,i2)) );
            end
        end
    end
    % fisher transformation
    z(z > 1-1e-15) = 1-1e-15; % to prevent this: atanh(1) = Inf
    z = atanh(z);

    % create output dir if necassary
    imgFname = fullfile(outDir,'vmhc',strcat(subjectName,'.nii'));
    [pth,~,~] = fileparts(imgFname);
    if ~exist(pth,'dir'), mkdir(pth); end
    % save fcunctional connectivity map
    afxVolumeWrite(imgFname,z,dim,'int16',mat,'vmhc: fisher transformed corrcoeff');
    fprintf('done\n')
end
