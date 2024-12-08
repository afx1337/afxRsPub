function f = afxPlotSections(fnamesNiiData,z,titles,threshData,limitsData,resample,dimension)
    dimCaption = {'x','y','z'};
    if nargin < 6, resample = 0; end
    if nargin < 7, dimension = 3; end

    %fnamesNiiTemplate = fullfile('templates','sch2better.nii');
    %limitsTemplate = [45 120];
    fnamesNiiTemplate = fullfile('templates','MNI152_T1_0.5mm_masked.nii');
    limitsTemplate = [0 255];

    % load template
    [Ytemplate,XYZmm] = spm_read_vols(spm_vol(fnamesNiiTemplate));
    Ytemplate = scaleDat(Ytemplate,limitsTemplate);
    dim = size(Ytemplate);
    if dimension == 2,     Ytemplate = permute(Ytemplate,[1 3 2]);
    elseif dimension == 1, Ytemplate = permute(Ytemplate,[2 3 1]);
    end
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    tmp = reshape(XYZmm(dimension,:),dim);
    if dimension == 3,     zMNI = squeeze(tmp(1,1,:));
    elseif dimension == 2, zMNI = squeeze(tmp(1,:,1));
    elseif dimension == 1, zMNI = squeeze(tmp(:,1,1));
    end
    clear tmp;

    % load data
    if ~iscell(fnamesNiiData)
        Ydata = getData(fnamesNiiData,XYZmm,dim,resample);
        if dimension == 2,     Ydata = permute(Ydata,[1 3 2]);
        elseif dimension == 1, Ydata = permute(Ydata,[2 3 1]);
        end
    end
        
    nZ = length(z);
    %n = max(round(sqrt(nZ)*1.5));
    %m = ceil(nZ/n);
    m = 1; n = nZ;

    f = figure;
    for iZ = 1:length(z)
        if iscell(fnamesNiiData)
            % load data
            Ydata = getData(fnamesNiiData{iZ},XYZmm,dim,resample);
            if dimension == 2,     Ydata = permute(Ydata,[1 3 2]);
            elseif dimension == 1, Ydata = permute(Ydata,[2 3 1]);
            end
            curTitle = titles{iZ};
        else
            curTitle = sprintf('%s = %d',dimCaption{dimension},z(iZ));
        end
        
        % get slice number
        if isnan(z(iZ)) % maximum value
            tmp = squeeze(max(max(Ydata)));
            curZMNI = find(tmp == max(tmp),1); clear tmp;
        else % coordinate
            curZMNI =  find(abs(zMNI -  z(iZ)) - min(abs(zMNI -  z(iZ))) == 0);
        end
        
        % plot image
        subplot(m,n,iZ);
        plotNii(Ytemplate,Ydata,threshData,curZMNI,limitsData,dimension,z(iZ));
        title(curTitle)
    end

end

function plotNii(Ytemplate,Ydata,threshData,z,limitsData,dimension,mni)
    % template slice
    Ytemplate = afxUpsample2d(Ytemplate(:,end:-1:1,z)');
    % white background
%     tmp = bwconncomp(Ytemplate < .1);
%     for i = 1:length(tmp.PixelIdxList)
%         if nnz(tmp.PixelIdxList{i} == 1) > 0
%             Ytemplate(tmp.PixelIdxList{i}) = 1;
%             break;
%         end
%     end
    Ytemplate = afxwWhiteBG(Ytemplate);
        
    % data slice
    %Ydata = afxUpsample2d(Ydata(:,end:-1:1,z)');
    Ydata = afxUpsample2d(Ydata(:,end:-1:1,z)');
    Ydata(Ydata < threshData(1)) = NaN;
    Ydata(Ydata > threshData(2)) = NaN;
    mask = isnan(Ydata)| Ydata == 0;
    
    if limitsData(1) == -Inf
        limitsData(1) = min(Ydata(:));
    elseif limitsData(2) == Inf
        limitsData(2) = max(Ydata(:));
    end
    Ydata = scaleDat(Ydata,limitsData);
    
    % mask
    Ydata(mask) = 0;
    Ytemplate(~mask) = 0;
    
    % data to rgb
    if threshData(2) < 0
        cMap = winter(256);
        cMap = cMap(end:-1:1,:);
    else
        cMap = hot(256);
    end
    img = grs2rgb(Ydata,cMap(32:end,:));
    
    % rgb data + overlay
    img(:,:,1) = img(:,:,1).*~mask+Ytemplate;
    img(:,:,2) = img(:,:,2).*~mask+Ytemplate;
    img(:,:,3) = img(:,:,3).*~mask+Ytemplate;

    % slight smoothing (antialiasing)
    img(:,:,1) = afxSmooth2d(img(:,:,1));
    img(:,:,2) = afxSmooth2d(img(:,:,2));
    img(:,:,3) = afxSmooth2d(img(:,:,3));

    if dimension == 1 && mni <= 0
        img = img(:,end:-1:1,:);
    end
    
    % show image
    imshow(img);
end

function img = afxSmooth2d(img)
 %[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 %f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 %f=f./sum(f(:));

    a = 1;
    b = .25;
    c = .062;
    kernel = [c b c; b a b; c b c];
    kernel = kernel/sum(kernel(:));
    img(2:end-1,2:end-1) = 0 ...
        + img(2:end-1,2:end-1)*kernel(5) ...
        + img(1:end-2,2:end-1)*kernel(2) ...
        + img(2:end-1,1:end-2)*kernel(2) ...
        + img(3:end-0,2:end-1)*kernel(2) ...
        + img(2:end-1,3:end-0)*kernel(2) ...
        + img(1:end-2,1:end-2)*kernel(1) ...
        + img(1:end-2,3:end-0)*kernel(1) ...
        + img(3:end-0,1:end-2)*kernel(1) ...
        + img(3:end-0,3:end-0)*kernel(1);
end

function dat = afxwWhiteBG(dat)
    datBg = dat < 0.1;
    bgLabels = zeros(size(datBg));
    bgLabels(1,:) = 1; bgLabels(end,:) = 1; bgLabels(:,1) = 1; bgLabels(:,end) = 1;
    dim = size(dat);
    cnt = 1;
    for x = 2:dim(2)-1
        for y = 2:dim(1)-1
            coordx = [ x     x-1 ];
            coordy = [ y-1   y   ];
            tmp1 = [datBg(coordy(1),coordx(1))    datBg(coordy(2),coordx(2))];
            tmp2 = [bgLabels(coordy(1),coordx(1)) bgLabels(coordy(2),coordx(2))];
            curLabels = tmp2(tmp1);
            curLabels = curLabels(curLabels>0);
            if isempty(curLabels)
                bgLabels(y,x) = cnt;
                cnt = cnt + 1;
            else
                bgLabels(y,x) = min(curLabels);
                if length(curLabels) > 1 && curLabels(1) ~= curLabels(2)
                    bgLabels(bgLabels == max(curLabels)) = min(curLabels);
                end
            end
        end
    end
    dat(bgLabels == 1) = 1;
    dat(1,:) = 1;
    dat(end,:) = 1;
    dat(:,1) = 1;
    dat(:,end) = 1;
end

function dat = scaleDat(dat,limits)
    if isempty(limits)
        limits(1) = min(dat(:));
        limits(2) = max(dat(:));
    end
    dat = (dat - limits(1))./diff(limits);
    dat(dat < 0) = 0;
    dat(dat > 1) = 1;
end

function dat = getData(fname,XYZ,dim,resample)
	dat = afxVolumeResample(fname,XYZ,resample);
    dat = reshape(dat,dim);
end

function img = afxUpsample2d(img)
    % resample data to higher resolution
    %img = imresize(double(img),2,'bilinear');
    %img = afxUpsample1d(afxUpsample1d(img')');
    img = afxUpsample1d(afxUpsample1d(img')');
end

function Y = afxUpsample1d(X)
    Y = nan(size(X,1)*2-1,size(X,2));
    for iy = 1:size(X,1)-1
        Y(2*iy-1,:) = X(iy,:);
        Y(2*iy,:) = (X(iy,:)+X(iy+1,:))/2;
    end
    Y(end,:) = X(end,:);
end


function res = grs2rgb(img, map)
    %%Convert grayscale images to RGB using specified colormap.
    %	IMG is the grayscale image. Must be specified as a name of the image
    %	including the directory, or the matrix.
    %	MAP is the M-by-3 matrix of colors.
    %
    %	RES = GRS2RGB(IMG) produces the RGB image RES from the grayscale image IMG
    %	using the colormap HOT with 64 colors.
    %
    %	RES = GRS2RGB(IMG,MAP) produces the RGB image RES from the grayscale image
    %	IMG using the colormap matrix MAP. MAP must contain 3 columns for Red,
    %	Green, and Blue components.
    %
    %	Example 1:
    %	open 'image.tif';
    %	res = grs2rgb(image);
    %
    %	Example 2:
    %	cmap = colormap(summer);
    % 	res = grs2rgb('image.tif',cmap);
    %
    % 	See also COLORMAP, HOT
    %
    %	Written by
    %	Valeriy R. Korostyshevskiy, PhD
    %	Georgetown University Medical Center
    %	Washington, D.C.
    %	December 2006
    %
    % 	vrk@georgetown.edu

    % Check the arguments
    if nargin<1
        error('grs2rgb:missingImage','Specify the name or the matrix of the image');
    end;

    if ~exist('map','var') || isempty(map)
        map = hot(64);
    end;

    [l,w] = size(map);

    if w~=3
        error('grs2rgb:wrongColormap','Colormap matrix must contain 3 columns');
    end;

    if ischar(img)
        a = imread(img);
    elseif isnumeric(img)
        a = img;
    else
        error('grs2rgb:wrongImageFormat','Image format: must be name or matrix');
    end;

    % Calculate the indices of the colormap matrix
    a = double(a);
    a(a==0) = 1; % Needed to produce nonzero index of the colormap matrix
    ci = ceil(l*a/max(a(:)));

    % Colors in the new image
    [il,iw] = size(a);
    r = zeros(il,iw);
    g = zeros(il,iw);
    b = zeros(il,iw);
    r(:) = map(ci,1);
    g(:) = map(ci,2);
    b(:) = map(ci,3);

    % New image
    res = zeros(il,iw,3);
    res(:,:,1) = r;
    res(:,:,2) = g;
    res(:,:,3) = b;
end