function alignRetinaPatch(bigPath, smallPath, varargin)
% alignRetinaPatch V2 —— manual / read existing tform，align high-resolution small sub-region images with labelled whole retina image；
%                     and only export the overlapped area
% e.g. alignRetinaPatch('retina_full.tif', 'patch.tif', ...
%                'LoadTform','patch_tform.mat', ...
%                'Blend',true, 'Alpha',0.4, ...
%                'OverlapOnly',true, ...
%                'SaveWarped','patch_with_grid.tif');
% main inputs
% ----
%   bigPath   : path for whole retina image labelled with RGC positions
%   smallPath : path for smaller sub-regions with clear visible cells
%   matPath   : path for the neurons' coordinates in pixel value for the
%   whole retina image
%
% optional inputs
% -----------------
%   'LoadTform'   : ''(d) | 'xxx.mat'  if exists load tform_total / tform
%   'OverlapOnly' : false(d) | true    only export the overlapped area
%
% 'TransformType': default: projective(project the feature with the ),'AutoRefine','SaveWarped',
% 'Blend','Alpha','SaveComposite','CropRect'……

% ---------- 0. setting parameters ----------
p = inputParser;
addParameter(p,'TransformType','projective');
addParameter(p,'AutoRefine',false);
addParameter(p,'SaveWarped','');
addParameter(p,'Blend',true);
addParameter(p,'Alpha',0.5);
addParameter(p,'LoadTform','');        
parse(p,varargin{:});

tfmType   = p.Results.TransformType;
doRefine  = p.Results.AutoRefine;
saveWarp  = p.Results.SaveWarped;
doBlend   = p.Results.Blend;
alphaWgt  = p.Results.Alpha;
loadTform = p.Results.LoadTform;

if isempty(saveWarp)
    saveWarp = [extractBefore(smallPath,'.tif') '_matched.tif'];
end
if isempty(loadTform)
    loadTform = [extractBefore(smallPath,'.tif') '_tform.mat'];
end
% ---------- 1. read images ----------
big   = imread(bigPath);
small = imread(smallPath);
RS = imref2d(size(small));

% ---------- 2. extract tform ----------
if ~isempty(loadTform) && isfile(loadTform)
    % ---- 2‑A  direct import from saved .mat file ----
    S = load(loadTform);
    if isfield(S,'tform_total')
        tform_total = S.tform_total;
    elseif isfield(S,'tform')
        tform_total = S.tform;
    else
        error('No variable tform_total / tform found in %s',loadTform);
    end
    fprintf('Loaded existing transform from %s\n',loadTform);
else
    % ---- 2‑B  manual alignment ----
    figure('Name','Choose Correspondences','NumberTitle','off');
    subplot(1,2,1); imshow(big,[]);   title('BIG — pick 4+');
    subplot(1,2,2); imshow(small,[]); title('SMALL — pick same 4+');
    fprintf(['\nStep 1  left➜right points alignment，≥4 pairs；Press ENTER when finished \n']);
    bigPts = []; smallPts = [];
    while true
        subplot(1,2,1); h1 = drawpoint('Color','y','MarkerSize',8);
        subplot(1,2,2); h2 = drawpoint('Color','y','MarkerSize',8);
        bigPts   = [bigPts;   h1.Position];
        smallPts = [smallPts; h2.Position];
        if size(bigPts,1) >= 4
            cont = questdlg('Need more points?','Continue?','Yes','No','No');
            if strcmp(cont,'No'), break; end
        end
    end
    close;
    % initialise tform
    tform = fitgeotrans(smallPts, bigPts, tfmType);
    RB = imref2d(size(big));
    smallWarp = imwarp(small, tform,'OutputView',RB);
    % automatic refine（optional）
    if doRefine
        [opt,met] = imregconfig('multimodal');
        tInc = imregtform(smallWarp,big,'affine',opt,met);
        tform_total = affine2d(tInc.T * tform.T);
    else
        tform_total = tform;
    end
    % save tform_total for next time if satisfied
    save([extractBefore(smallPath,'.tif') '_tform.mat'],'tform_total');
    fprintf('Transform saved as tform_total.mat\n');
end

RB = imref2d(size(big));
small2big = imwarp(small, tform_total, 'OutputView', RB);

figure('Name','Overlay Check (big‑space)');        % for preview
imshowpair(big, small2big,'blend');
title('Close window to continue'); uiwait(gcf); close;

%% warp RGB pic→ PATCH coords
bigRGBsmall = imwarp(big, invert(tform_total), 'OutputView', RS);

%% extract the grey base image（after warp greyscale from big image）+ smaller image α blend
bigGraySmall = rgb2gray(bigRGBsmall);

maskPatch = small > 0;
classType = class(small);

blend_d = double(bigGraySmall);
if doBlend
    blend_d(maskPatch) = (1-alphaWgt)*double(bigGraySmall(maskPatch)) ...
                       + alphaWgt        *double(small(maskPatch));
else
    blend_d(maskPatch) = double(small(maskPatch));
end
baseGray = cast(blend_d,classType);               % uint16 grey scale

%% Find color mask in the big pic
R = bigRGBsmall(:,:,1); G = bigRGBsmall(:,:,2); B = bigRGBsmall(:,:,3);
colorMaskSmall = abs(int16(R)-int16(G))>5 | ...
                 abs(int16(R)-int16(B))>5 | ...
                 abs(int16(G)-int16(B))>5;    
colorMaskSmall = bwareaopen(colorMaskSmall,20);   

%% combine 
% 
base8  = uint8(255*mat2gray(baseGray));
compRGB = repmat(base8,1,1,3);

% 
for ch = 1:3
    tmp            = compRGB(:,:,ch);
    A = bigRGBsmall(:,:,ch);
    tmp(colorMaskSmall) = A(colorMaskSmall);
    compRGB(:,:,ch) = tmp;
end

%% only keep small image OR colored outline
keepMask = maskPatch | colorMaskSmall;
[r,c] = find(keepMask);
compRGB = compRGB(min(r):max(r), min(c):max(c), :);

%% save
imwrite(compRGB, saveWarp,'Compression','none');
fprintf('Saved patch with ORIGINAL COLOR frames → %s\n',saveWarp);

end
