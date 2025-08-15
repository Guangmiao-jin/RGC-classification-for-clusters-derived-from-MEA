function chipGridAlignspikedensity(imgPath, matPath,varargin)

p = inputParser;
addParameter(p, 'BinWidth', 20, @isnumeric); % 
addParameter(p, 'SmoothWindow', 3, @isnumeric); % 
parse(p, varargin{:});

%% read the tiff file
I = imread(imgPath);
[H,W] = size(I);
figure('Name','Chip Grid Calibration','NumberTitle','off');
hAx = axes; imshow(I,[],'Parent',hAx); title('Step 1: Zoom & pick the FOUR corner channels, top left -> top right -> bottom left -> bottom right');

%% 1. extract the centres of four corners
cornerPt = [];                 % 4×2
for k = 1:4
    h = drawpoint('Color','r','MarkerSize',10,'Label',sprintf('%d',k));
    cornerPt(k,:) = h.Position;
end
%% 2. calculate the grids' distances
% assume cornerPt goes in the following order：c00, c0N, cM0, cMN
c00 = cornerPt(1,:);         % top left
c0N = cornerPt(2,:);         % top right
cM0 = cornerPt(3,:);         % bottom left

% calculate the distance between each channels horizontally and vertically
ux = (c0N - c00)/(64-1);     % horizontal distance
uy = (cM0 - c00)/(64-1);     % vertical distance

[Xidx,Yidx] = meshgrid(0:63,0:63); 
gridX = c00(1) + ux(1)*Xidx + uy(1)*Yidx;
gridY = c00(2) + ux(2)*Xidx + uy(2)*Yidx;

%% 3. visualization
hold on;
scatter(gridX(:),gridY(:),15,'y','filled');   % grids' points
% draw the borders
plot([gridX(1,1) gridX(1,end)],[gridY(1,1) gridY(1,end)],'y-','LineWidth',1.5);
plot([gridX(end,1) gridX(end,end)],[gridY(end,1) gridY(end,end)],'y-','LineWidth',1.5);
plot([gridX(1,1) gridX(end,1)],[gridY(1,1) gridY(end,1)],'y-','LineWidth',1.5);
plot([gridX(1,end) gridX(end,end)],[gridY(1,end) gridY(end,end)],'y-','LineWidth',1.5);
title('Completed！Yellow dosts indicate 64×64 channels'); pause;
close(gcf);

%% 4. coordinates output
coords = [reshape(Xidx,[],1) reshape(Yidx,[],1) reshape(gridX,[],1) reshape(gridY,[],1)];
assignin('base','channelCoords',coords);  % saved to workspace
disp('Variable channelCoords = [row col x y] has been saved to the workspace');

%% 5. read _cluster.mat file and annotate spike density on tiff image
S = load(matPath);          
col_f = round(S.data.spikeX/60);    % 0-63
row_f = round(S.data.spikeY/60);    % 0-63
valid = col_f>=0 & col_f<=64 & row_f>=0 & row_f<=64;
col_f = col_f(valid); 
row_f = row_f(valid);

Xp = [];
Yp = [];
k = find(valid==1);
for x = 1:length(k)
    Xp = [Xp; gridX(64-row_f(k(x)),col_f(k(x))+1)];   % → pixel x
    Yp = [Yp; gridY(64-row_f(k(x)),col_f(k(x))+1)];  % → pixel y
end
%% 2D histogram 
cMN = c00 + 63*ux + 63*uy;  
xMin = max(1, floor(min([c00(1), c0N(1), cM0(1), cMN(1)])));
xMax = min(W, ceil( max([c00(1), c0N(1), cM0(1), cMN(1)]) ));
yMin = max(1, floor(min([c00(2), c0N(2), cM0(2), cMN(2)])));
yMax = min(H, ceil( max([c00(2), c0N(2), cM0(2), cMN(2)]) ));

Xedges = xMin:(xMax-xMin)/64:xMax; if Xedges(end)<xMax, Xedges=[Xedges xMax]; end
Yedges = yMin:(yMax-yMin)/64:yMax; if Yedges(end)<yMax, Yedges=[Yedges yMax]; end

[countsXY,~,~] = histcounts2(Xp, Yp, Xedges, Yedges);   % [Nx × Ny]
% (row=Y，column=X）
C = countsXY';  % Ny × Nx

% generate meshgrid of bin centres
xCenters = (Xedges(1:end-1) + Xedges(2:end))/2;
yCenters = (Yedges(1:end-1) + Yedges(2:end))/2;
[XcBins, YcBins] = meshgrid(xCenters, yCenters);  

% four-corner mask
polyX = [c00(1) c0N(1) cMN(1) cM0(1)];
polyY = [c00(2) c0N(2) cMN(2) cM0(2)];
inChip = inpolygon(XcBins, YcBins, polyX, polyY);
C(~inChip) = 0;

% —— log10（no spikes-> transparent）——
nz = C(C>0);
if isempty(nz)
    lo = 1; hi = 1;          
else
    lo = max(1, prctile(nz, 1));     % 1% min
    hi = prctile(nz, 100);            % 100% max
end
Clog = zeros(size(C), 'like', C);
maskNZ = C>0;
Clog(maskNZ) = log10(C(maskNZ));
CLim = log10([lo hi]);

% transparency
Alpha = zeros(size(C));
Alpha(maskNZ) = 1;  % 0..1
Alpha = min(max(Alpha,0),1);

%% 6. plotting
I8 = im2uint8(rescale(I));
Irgb = repmat(I8, 1, 1, 3);   % H×W×3 grey RGB
figure('Name','Chip + 2D histogram (bin=50px, log scale)','NumberTitle','off');
image([1 W],[1 H], Irgb); axis image; hold on;

hImg = imagesc(xCenters, yCenters, Clog);    % log10 values
set(hImg, 'AlphaData', Alpha);

n = 256;  % colormap 分辨率
baseMap = parula(512);   % 或 turbo(512)，分辨率更高便于裁剪

% parula: 前段蓝->中段绿->黄；这里取前 2/3 区间
idxRange = round(linspace(1, round(size(baseMap,1)*2/3), n));
cmap = baseMap(idxRange, :);

colormap(gca, cmap);
clim(CLim);  % 或 caxis(CLim)


% colormap convert to RGB 
cmapN = size(cmap,1);

normVal = (Clog - CLim(1)) / (CLim(2) - CLim(1));
normVal(~maskNZ) = NaN;  
idxVal = round(normVal * (cmapN-1) + 1);
idxVal = max(1, min(cmapN, idxVal)); 

[Xpix, Ypix] = meshgrid(1:W, 1:H);
overlayR = nan(H, W);
overlayG = nan(H, W);
overlayB = nan(H, W);
alphaPix = nan(H, W);

overlayR(:) = interp2(xCenters, yCenters, reshape(cmap(idxVal,1),size(C)), Xpix, Ypix, 'nearest');
overlayG(:) = interp2(xCenters, yCenters, reshape(cmap(idxVal,2),size(C)), Xpix, Ypix, 'nearest');
overlayB(:) = interp2(xCenters, yCenters, reshape(cmap(idxVal,3),size(C)), Xpix, Ypix, 'nearest');
alphaPix(:) = interp2(xCenters, yCenters, Alpha, Xpix, Ypix, 'nearest');

overlayR(isnan(overlayR)) = 0;
overlayG(isnan(overlayG)) = 0;
overlayB(isnan(overlayB)) = 0;
alphaPix(isnan(alphaPix)) = 0;

IrgbOut = uint8( ...
    (1 - alphaPix).*double(Irgb) + ...
    alphaPix.*cat(3, overlayR*255, overlayG*255, overlayB*255) );

close all;

figure('Name','Full-res Spike Density Overlay','NumberTitle','off');
imshow(IrgbOut); axis image; hold on;
colormap(gca, cmap);
clim(CLim);  
cb = colorbar;

tickVals = [];
for k = floor(log10(lo)) : ceil(log10(hi))
    tickVals = [tickVals, 10^k]; 
end
tickVals = tickVals(tickVals>=lo & tickVals<=hi);
cb.Ticks = log10(tickVals);
cb.TickLabels = string(log10(tickVals));
ylabel(cb, 'log10 (Spike count)');

title('Spike density graph');

%% 7. saving graph
outPath_fig =[extractBefore(imgPath,'.tif') '_spikedensity_colormap.tif'];
exportgraphics(gcf, outPath_fig, 'Resolution', 300);  % high-res graph
fprintf('✓ Figure with colorbar saved → %s\n', outPath_fig);

close all;

end

