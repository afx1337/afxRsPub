function afxMirror(niiIn)
    % afxMirror(niiIn)
    % max.wawrzyniak@medizin.uniklinik-leipzig.de

    if nargin < 1, niiIn = spm_select([1 Inf],'image'); end

	% default niiFlip

    
    for i = 1:size(niiIn,1)
        curNii = niiIn(i,:);
        [path, name, ext] = fileparts(curNii);
        niiFlip = fullfile(path,strcat(name,'-mirrored',ext(1:4)));
    
        % load nii
        V = nifti(curNii);
        % save data
        dat = V.dat(:,:,:);
        % create new empty nifti
        V.dat.fname = niiFlip;
        V.dat(:,:,:) = 0;
        create(V);
    
        % write flipped data
        V.dat(:,:,:) = dat(V.dat.dim(1):-1:1,:,:);
    end
end