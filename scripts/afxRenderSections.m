function montageRGB = afxRenderSections(templateVol, templateMat, sliceSpec, sliceList, overlays, opts, outFilename)
% afxRenderSections  Render multiple MNI slices and concatenate them into a montage.
%
%   montageRGB = afxRenderSections(templateVol, templateXYZmm, sliceSpec, sliceList, overlays, opts [, outFilename])
%
%   Renders a list of slices from a 3D template volume in MNI space and
%   applies optional overlays. Each slice is rendered by renderMNISlice,
%   and the resulting images are concatenated into a montage.
%
%   The function returns an RGB montage image (double, range [0 1]).
%
% -------------------------------------------------------------------------
% Inputs
% -------------------------------------------------------------------------
%   templateVol (3D numeric array)
%       Template volume (e.g., MNI normalized anatomical image).
%
%   templateXYZmm (3 x N numeric array)
%       Voxel-to-world mapping (MNI coordinates). N = numel(templateVol).
%       Rows correspond to x, y, z coordinates.
%
%   sliceSpec (struct)
%       Slice specification. Must contain:
%           .orientation       - 'axial', 'coronal', or 'sagittal'
%       sliceSpec.mm is overwritten for each slice in sliceList.
%
%   sliceList (vector)
%       List of slice positions in mm along the selected orientation axis.
%       Example (axial): [-30 -10 10 30]
%
%   overlays (struct array)
%       Overlay definitions. Each element must contain:
%           .vol               - overlay volume
%           .mat               - overlay affine mat
%           .fname             - overlay file name (if vol/mat is empty)
%           .thr               - [low high] threshold (outside set to 0)
%           .mode              - 'solid' | 'outline' | 'outline+solid' | 'gradient'
%           .color             - [r g b] (0..1) for solid/outline
%           .alpha             - transparency (0..1) for solid/gradient
%           .relativeLineWidth - outline thickness relative to diagonal
%           .colormap          - string or Nx3 colormap for gradient
%           .clim              - [low high] for gradient mapping
%
%   opts (struct, optional)
%       Rendering and montage options. Fields (defaults shown):
%           .intensityWindow  - [low high] linear intensity window for underlay
%           .scaleFactor      - 1          Upsampling factor for slices
%           .interp           - 'bilinear'/'nearest' 
%                                          Interpolation for upsampling
%           .maskingThrBg     - 0.15       Threshold for masking background
%           .maskOverlays     - true       Apply mask to overlays as well
%           .showLabel        - true       Add label row above each slice
%           .labelFontSize    - 14         Font size for labels
%           .labelColor       - 'black'    Text color
%           .labelPosition    - 'top'      'top' or 'bottom' of slice block
%           .nCols            - []         Number of montage columns.
%                                          If empty, a 16:9-like layout is used.
%
%   outFilename (optional)    - 'img.png'  Render to file directly
%
% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
%   montageRGB (HxWx3 numeric array)
%       RGB montage image (double, range [0 1]). Each slice is rendered
%       with overlays and combined with a label row (if showLabel=true).
%
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
%   - templateXYZmm is used to compute the MNI coordinate axis only once
%     for performance (not per slice).
%   - renderMNISlice performs the actual rendering of a single slice.
%   - The montage is assembled row-major (left-to-right, top-to-bottom).
%
% -------------------------------------------------------------------------
% Example
% -------------------------------------------------------------------------
% [templateVol,XYZmm,dim] = afxLoadFunc(char('underlay.nii'));
% templateVol = reshape(templateVol, dim);
% dat = afxVolumeResample('overlay.nii',XYZmm,1);
% dat = reshape(dat,dim);
% 
% sliceSpec.orientation = 'axial';
% sliceList = -40:20:60;
% 
% opts.intensityWindow = [0 255];
% opts.interp = 'bilinear';
% opts.scaleFactor = 2;
% 
% overlay(1).vol = dat;
% overlay(1).thr = [-Inf Inf];
% overlay(1).mode = 'gradient';
% overlay(1).colormap = 'hot';
% overlay(1).alpha = .6;
% overlay(1).clim = [-3 3];
% 
% montageRGB = afxRenderSections(templateVol, XYZmm, sliceSpec, sliceList, overlay, opts);
% imshow(montageRGB);

if nargin < 6, opts = struct(); end
if ~isfield(opts,'scaleFactor'), opts.scaleFactor = 1; end
if ~isfield(opts,'interp'), opts.interp = 'bilinear'; end
if ~isfield(opts,'maskingThrBg'), opts.maskingThrBg = 0.15; end
if ~isfield(opts,'maskOverlays'), opts.maskOverlays = true; end
if ~isfield(opts,'nCols'), opts.nCols = []; end
if ~isfield(opts,'showLabel'), opts.showLabel = true; end
if ~isfield(opts,'labelFontSize'), opts.labelFontSize = 14; end
if ~isfield(opts,'labelPosition'), opts.labelPosition = 'top'; end
if ~isfield(opts,'labelColor'), opts.labelColor = 'black'; end

opts.labelFontSize = opts.labelFontSize * opts.scaleFactor;

n = numel(sliceList);
blocks = cell(n,1);

axStr = {'x','y','z'};

MNImm = prepareMNICoords(size(templateVol), templateMat, sliceSpec);

% load overlay data if necessary 
for i = 1:numel(overlays)
   if isfield(overlays,'fname') && ~isempty(overlays(i).fname) && (~isfield(overlays,'vol') || isempty(overlays(i).vol))
       %overlays(i).vol = reshape(afxVolumeResample(overlays(i).fname,templateXYZmm,0),size(templateVol));
       [ovVol,~,ovDim,ovMat] = afxLoadFunc(char(overlays(i).fname));
       if min(ovVol(:)) == max(ovVol(:)), ovVol(:) = 0; end
       ovVol(isnan(ovVol(:))) = 0;
       overlays(i).vol = reshape(ovVol,ovDim);
       overlays(i).mat = ovMat;
   end
end

for i = 1:n
    sliceSpec.mm = sliceList(i);
    
    % render slice without label
    [rgbSlice, info] = renderMNISlice(templateVol, templateMat, MNImm, sliceSpec, overlays, opts);

    % create label string
    if opts.showLabel
        labelStr = sprintf('%s = %d', axStr{info.sliceAxis}, info.mm);
        % label row
        [~, w, ~] = size(rgbSlice);
        labelRow = renderLabelRow(labelStr, w, opts);
    else
        labelRow = [];
    end

    % combine
    if strcmp(opts.labelPosition,'bottom')
        blocks{i} = [rgbSlice; labelRow];
    else
        blocks{i} = [labelRow; rgbSlice];
    end
end

montageRGB = concatImages(blocks, opts);

if exist('outFilename','var')
    [pth,~,~] = fileparts(outFilename);
    if ~exist(pth,'dir'), mkdir(pth); end
    imwrite(montageRGB,outFilename);
end

end


function montageRGB = concatImages(imgs, opts)

n = numel(imgs);
[h,w,~] = size(imgs{1});
imgRatio = w / h;  % Seitenverh�ltnis des ersten Bildes

% determine layout
if isempty(opts.nCols)
    targetRatio = 16/9 / imgRatio;
    nCols = ceil(sqrt(n * targetRatio));
    opts.nCols = nCols;
end
nCols = opts.nCols;
nRows = ceil(n / nCols);

montageRGB = ones(nRows*h, nCols*w, 3);

for k = 1:n
    r = floor((k-1)/nCols);
    c = mod((k-1), nCols);

    montageRGB(r*h+(1:h), c*w+(1:w), :) = imgs{k};
end

end


function MNImm = prepareMNICoords(dim, mat, sliceSpec)

orientation = lower(sliceSpec.orientation);

% ---------------------------
% Which world axis is sliced?
% ---------------------------
switch orientation
    case 'axial'
        sliceWorldAxis = 3; % Z
    case 'coronal'
        sliceWorldAxis = 2; % Y
    case 'sagittal'
        sliceWorldAxis = 1; % X
    otherwise
        error('Unknown orientation: %s', orientation);
end

% ---------------------------
% MNI 
% ---------------------------

MNImm = getSliceCoordinates(mat, dim, sliceWorldAxis);

end


function [rgb, info] = renderMNISlice(templateVol, matWorld, MNImm, sliceSpec, overlays, opts)
% ---------------------------
% Base slice
% ---------------------------
sliceSpec.interp = 'bilinear';
sliceSpec.scaleFactor = opts.scaleFactor;
[sliceT, info] = extractVoxelSlice(templateVol, matWorld, size(templateVol), matWorld, MNImm, sliceSpec);
[rgb,bgMask] = renderBaseSlice(sliceT, opts.intensityWindow, opts.maskingThrBg);
sliceSpec.interp = opts.interp;

if ~opts.maskOverlays
    bgMask(:) = 0;
end

% ---------------------------
% Overlays
% ---------------------------
for i = 1:numel(overlays)
    ov = overlays(i);

    % defaults
    if ~isfield(ov,'color') || isempty(ov.color), ov.color = [.7 0 0]; end
    if ~isfield(ov,'alpha') || isempty(ov.alpha), ov.alpha = 1; end
    if ~isfield(ov,'relativeLineWidth') || isempty(ov.relativeLineWidth), ov.relativeLineWidth = 0.002; end
    if ~isfield(ov,'colormap') || isempty(ov.colormap), ov.colormap = 'hot'; end
    if ~isfield(ov,'thr') || isempty(ov.thr), ov.thr = [-Inf Inf]; end
    if ~isfield(ov,'clim') || isempty(ov.clim)
        tmp = ov.vol(ov.vol(:) ~= 0);
        ov.clim = [max(min(tmp), ov.thr(1)) min(max(tmp), ov.thr(2))];
    end
    
    % Extract overlay slice
    sliceO = extractVoxelSlice(ov.vol, ov.mat, size(templateVol), matWorld, MNImm, sliceSpec);
    sliceO(sliceO < ov.thr(1) | sliceO > ov.thr(2)) = 0;

    switch lower(ov.mode)
        case 'solid'
            rgb = applyOverlaySolid(rgb, bgMask, sliceO, struct( ...
                'color', ov.color, 'alpha', ov.alpha));

        case 'outline'
            rgb = applyOverlayOutline(rgb, bgMask, sliceO, struct( ...
                'color', ov.color, 'relativeLineWidth', ov.relativeLineWidth));

        case 'outline+solid'
            rgb = applyOverlaySolid(rgb, bgMask, sliceO, struct( ...
                'color', ov.color, 'alpha', ov.alpha));
            rgb = applyOverlayOutline(rgb, bgMask, sliceO, struct( ...
                'color', ov.color, 'relativeLineWidth', ov.relativeLineWidth));

        case 'gradient'
            rgb = applyOverlayGradient(rgb, bgMask, sliceO, struct( ...
                'colormap', ov.colormap, 'clim', ov.clim, 'alpha', ov.alpha));

        otherwise
            error('Unknown overlay mode: %s', ov.mode);
    end
end

% optional Gauss smoothing of alpha mask (anti-aliasing)
rgb = imgaussfilt(rgb, .6);

end


function [rgb,bgMask] = renderBaseSlice(slice2d, intensityWindow, cutoff)
% renderBaseSlice
% Render a grayscale RGB image from a 2D slice using strict linear
% intensity windowing, and remove connected background regions.
%
% bgOpts (optional struct):
%   .cutoff            - intensity threshold for candidate background

low  = intensityWindow(1);
high = intensityWindow(2);

if low >= high
    error('Intensity window must satisfy low < high.');
end

% ---------------------------
% Clamp
% ---------------------------
slice2d(slice2d < low)  = low;
slice2d(slice2d > high) = high;

% ---------------------------
% Scale to [0 1]
% ---------------------------
slice2d = (slice2d - low) / (high - low);

% ---------------------------
% Handle NaNs explicitly
% ---------------------------
slice2d(isnan(slice2d)) = 0;

% ---------------------------
% Background removal via connected components
% ---------------------------
[h, w] = size(slice2d);

candMask = slice2d <= cutoff;   % candidate background pixels
cc = bwconncomp(candMask, 8);
stats = regionprops(cc, slice2d, {'Area','MeanIntensity','PixelIdxList'});

bgMask = false(h,w);

for i = 1:cc.NumObjects
    idx = stats(i).PixelIdxList;
    %meanVal = stats(i).MeanIntensity;
    %varVal  = var(double(slice2d(idx)));

    % border region?
    [r,c] = ind2sub([h,w], idx);
    if any(r==1 | r==h | c==1 | c==w)
        bgMask(idx) = true;
        continue;
    end
end

% Set background to white (not black!)
slice2d(bgMask) = 1;

% ---------------------------
% Expand to RGB
% ---------------------------
rgb = repmat(slice2d, [1 1 3]);
end


function rgb = applyOverlaySolid(rgb, bgMask, overlaySlice, opts)
% applyOverlaySolid
%   Apply a solid-color overlay to an RGB image.
%
%   overlaySlice must already be aligned to rgb (same orientation,
%   same transpose/flip).

mask = overlaySlice ~= 0 & ~bgMask;

if ~any(mask(:))
    return
end

alpha = opts.alpha;
color = reshape(opts.color, 1, 1, 3);

overlayRGB = repmat(color, size(rgb,1), size(rgb,2));

alphaMask = double(mask) * alpha;

rgb = blendOverlay(rgb, overlayRGB, alphaMask);
end


function rgb = applyOverlayOutline(rgb, bgMask, overlaySlice, opts)
% applyOverlayOutline
% Apply an outline overlay derived from a binary mask.
% Outline thickness is specified relative to the image diagonal.

mask = overlaySlice ~= 0 & ~bgMask;
if ~any(mask(:))
    return
end

% ---------------------------
% Perimeter
% ---------------------------
outline = bwperim(mask);

% ---------------------------
% Line width (relative)
% ---------------------------
[h, w] = size(mask);
diagPx = hypot(h, w);
lw = max(1, round(opts.relativeLineWidth * diagPx));

if lw > 1
    outline = imdilate(outline, strel('disk', lw));
end

% ---------------------------
% Build overlay RGB
% ---------------------------
overlayRGB = repmat(reshape(opts.color, 1, 1, 3), h, w);

% ---------------------------
% Alpha mask (opaque outline)
% ---------------------------
alphaMask = double(outline);

% ---------------------------
% Blend via zentrale Funktion
% ---------------------------
rgb = blendOverlay(rgb, overlayRGB, alphaMask);
end


function rgb = applyOverlayGradient(rgb, bgMask, overlaySlice, opts)
% applyOverlayGradient
% Apply a colormap-based gradient overlay with constant alpha.

clim = opts.clim;
if clim(1) >= clim(2)
    error('clim must satisfy low < high.');
end

% when 0 is not between data limits, perform implizit zero maskting
if ~(min(overlaySlice(:)) < 0 && max(overlaySlice(:)) > 0)
    bgMask = bgMask | overlaySlice == 0;
end
bgMask = bgMask | isnan(overlaySlice); % explicit nan mask

% ---------------------------
% Clamp & normalize
% ---------------------------
vals = overlaySlice;
vals(vals < clim(1)) = clim(1);
vals(vals > clim(2)) = clim(2);
vals = (vals - clim(1)) / (clim(2) - clim(1));

% ---------------------------
% Colormap
% ---------------------------
if ischar(opts.colormap) || isstring(opts.colormap)
    cmap = feval(opts.colormap, 256);
else
    cmap = opts.colormap;
end

if clim(2) < 0
    cmap = cmap(end:-1:1,:);
end

nC = size(cmap,1);
idx = round(vals * (nC-1)) + 1;
idx(bgMask) = 1; % dummy, won't be used

% ---------------------------
% Build overlay RGB
% ---------------------------
overlayRGB = reshape(cmap(idx,:), [size(overlaySlice,1), size(overlaySlice,2), 3]);

% ---------------------------
% Alpha mask (constant alpha where mask==true)
% ---------------------------
alphaMask = double(~bgMask) * opts.alpha;

% ---------------------------
% Blend via zentrale Funktion
% ---------------------------
rgb = blendOverlay(rgb, overlayRGB, alphaMask);

end


function rgbOut = blendOverlay(rgbBg, rgbOv, alphaMask)
% blendOverlay - blend overlay onto background using alpha mask
% alphaMask: 2D double [0..1]
% rgbOv, rgbBg: HxWx3 double [0..1]

% ensure range
alphaMask = max(0, min(1, alphaMask));

% make 3D
alpha3 = repmat(alphaMask, [1 1 3]);

% blend
rgbOut = alpha3 .* rgbOv + (1 - alpha3) .* rgbBg;
end


function labelRow = renderLabelRow(labelStr, w, opts)
% renderLabelRow - render label text into a row image (white background)

persistent fig ax

labelHeight = round(opts.labelFontSize * 2.5);
labelHeight = max(labelHeight, 24);

labelRow = ones(labelHeight, w, 3);

if isempty(labelStr)
    return;
end

% Figure einmal erstellen
if isempty(fig) || ~isvalid(fig)
    fig = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 w labelHeight]);
    ax = axes(fig, 'Position', [0 0 1 1], 'Units', 'pixels');
end

imshow(labelRow, 'InitialMagnification', 100);
axis off;

text(ax, w/2, labelHeight/2, labelStr, ...
    'Color', opts.labelColor, ...
    'FontSize', opts.labelFontSize, ...
    'Interpreter', 'none', ...
    'Units', 'pixels', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

frame = getframe(ax);
labelRow = im2double(frame.cdata);
close(fig);
end



function [slice2d, info] = extractVoxelSlice(vol, matVol, dimWorld, matWorld, MNImm, sliceSpec)

orientation = lower(sliceSpec.orientation);
mm          = sliceSpec.mm(:);

% ---------------------------
% Bounds check
% ---------------------------
if mm < min(MNImm) || mm > max(MNImm)
    error('Requested MNI slice (%.1f mm) outside volume.', mm);
end

% ---------------------------
% Which world axis is sliced?
% ---------------------------
switch orientation
    case 'axial'
        sliceWorldAxis = 3; % Z
    case 'coronal'
        sliceWorldAxis = 2; % Y
    case 'sagittal'
        sliceWorldAxis = 1; % X
    otherwise
        error('Unknown orientation: %s', orientation);
end

% ---------------------------
% Permute dimensions
% ---------------------------
switch orientation
    case 'axial'
        perm = [1 2 3];
    case 'coronal'
        perm = [1 3 2];
    case 'sagittal'
        perm = [2 3 1];
end
P = eye(4);
P(1:3,1:3) = P(1:3,perm);
vol    = permute(vol, perm);
matVol = matVol * P;
dimWorld = dimWorld(perm);
matWorld = matWorld * P;

[~, idx] = min(abs(MNImm - mm));

slice2d = resampleSliceAffine(dimWorld, matWorld, vol, matVol, idx, sliceSpec.interp);

slice2d = slice2d';
slice2d = slice2d(end:-1:1,:);

if sliceWorldAxis == 1 && mm <= 0
    slice2d = slice2d(:,end:-1:1);
end

[slice2d, ~] = upsampleSlice2D(slice2d, sliceSpec.scaleFactor, sliceSpec.interp);

% ---------------------------
% Info
% ---------------------------
info.voxelIndex  = idx;
info.mm          = MNImm(idx);
info.orientation = orientation;
info.sliceAxis   = sliceWorldAxis;
end

function [sliceUp, scaleInfo] = upsampleSlice2D(slice2d, scaleFactor, method)
% upsampleSlice2D - Upsample a 2D slice.
% scaleFactor: integer >= 1
% method: 'nearest' or 'bilinear'

if nargin < 3
    method = 'bilinear';
end

if scaleFactor <= 1
    sliceUp = slice2d;
    scaleInfo.scaleFactor = 1;
    return;
end

sliceUp = imresize(slice2d, scaleFactor, method);

scaleInfo.scaleFactor = scaleFactor;
end

function slice_resampled = resampleSliceAffine(dimWorld, matWorld, volLocal, matLocal, k0, interpMethod)
%RESAMPLESLICEAFFINE  Resample LowRes volume into HighRes slice space
%
% slice_resampled = resampleSliceAffine(Vh, Ah, Vl, Al, k0)
% slice_resampled = resampleSliceAffine(..., interpMethod)
%
% INPUT
%   Vh  : HighRes volume (Nx x Ny x Nz)
%   Ah  : 4x4 affine (High voxel -> MNI)
%   Vl  : LowRes volume  (Mx x My x Mz)
%   Al  : 4x4 affine (Low voxel -> MNI)
%   k0  : axial slice index in HighRes space
%   interpMethod : 'linear' (default) | 'nearest'
%
% OUTPUT
%   slice_resampled : LowRes volume sampled in HighRes slice space
%
% NOTE
%   MATLAB volume convention:
%   V(row, col, slice) = V(y, x, z)

if nargin < 6 || strcmp(interpMethod,'bilinear')
    interpMethod = 'linear';
end

% --- Size ---
Nhx = dimWorld(1);
Nhy = dimWorld(2);
Nhz = dimWorld(3);

if k0 < 1 || k0 > Nhz
    error('k0 outside HighRes volume range.');
end

% --- Transform High -> Low voxel space ---
T = inv(matLocal) * matWorld;

% --- HighRes slice grid ---
[i,j] = ndgrid(1:Nhx, 1:Nhy);
k = k0 * ones(size(i));

coords_high = [i(:)'; j(:)'; k(:)'; ones(1,numel(i))];

% --- Transform to LowRes voxel coordinates ---
coords_low = T * coords_high;

x_low = coords_low(1,:);
y_low = coords_low(2,:);
z_low = coords_low(3,:);

% --- Interpolator ---
F = griddedInterpolant(volLocal, interpMethod, 'none');

% IMPORTANT: MATLAB ordering (row=y, col=x, slice=z)
slice_resampled = F(x_low, y_low, z_low);

% reshape to 2D slice
slice_resampled = reshape(slice_resampled, Nhx, Nhy);

end


function z_coords = getSliceCoordinates(mat, volSize, sliceWorldAxis)

Nz = volSize(3);

% Wir nehmen die Voxelkoordinate im Zentrum der x/y Ebene
i0 = (volSize(1)+1)/2;
j0 = (volSize(2)+1)/2;

k = 1:Nz;

vox = [ ...
    i0*ones(1,Nz); ...
    j0*ones(1,Nz); ...
    k; ...
    ones(1,Nz)];

mni = mat * vox;

z_coords = mni(sliceWorldAxis,:);

end