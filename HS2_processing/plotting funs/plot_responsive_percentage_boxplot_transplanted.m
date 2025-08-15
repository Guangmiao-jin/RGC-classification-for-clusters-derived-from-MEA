function plot_responsive_percentage_boxplot_transplanted(baseFolder)
% plot_responsive_percentage_boxplot(baseFolder)
%
% go through all the dayxx subfiles in baseFolder
% read .mat file，calculate "responsive neuron %"，
% only keep those sample number ≥ 3 and plot boxplot。

%% 1. collect all data in each day file
groupFolders = dir(fullfile(baseFolder)); 
groupFolders = groupFolders([groupFolders.isdir] & ~ismember({groupFolders.name},{'.','..'}));
data      = [];      
counter    = 0;      % calculate how many days are counted

for i = 1:length(groupFolders)
    thisName = groupFolders(i).name;
    thisPath = fullfile(baseFolder, thisName);
    
    names   = groupFolders(i).name;

    % find .mat inside each file
    mats = dir(fullfile(thisPath, '*_totalneuronsV2.mat'));
    if numel(mats) < 3          % skip if less than 3 data points
        continue
    end

    % ----- extract the files -----
    respPerc = nan(numel(mats), 1);
    for j = 1:numel(mats)
        matObj = load(fullfile(thisPath, mats(j).name));

        % make sure the data is saved totalneurons structure
        if ~isfield(matObj, 'totalneurons')
            error('File %s lacks totalneurons structure', mats(j).name);
        end
        tn = matObj.totalneurons;

        % 1) find total neurons
        if ~isfield(tn, 'allneurons')
            error('File %s lacks allneurons ', mats(j).name);
        end
        totalN = numel(tn.allneurons);

        % 2) find unresponsive neurons
        notRespN = 0;
        if isfield(tn, 'num_notRespond')
            notRespN = numel(tn.num_notRespond);
        end

        % 3) calculate total percentages
        respPerc(j) = (totalN - notRespN) / totalN * 100;
    end

    % ----- save into data-----
    counter = counter + 1;
    data(counter).day           = names;
    data(counter).responsive    = respPerc;
    data(counter).nMice         = numel(respPerc);
end

%% 2. if there are no valid day，quit
if isempty(data)
    disp('⚠️  No day contains sample ≥ 3');
    return
end

%% 3. combine all the data，ready for plotting
allVals  = [];
allGroup = [];      
labels   = strings(numel(data),1);

for k = 1:numel(data)
    allVals            = [allVals;  data(k).responsive(:)];
    allGroup           = [allGroup; k * ones(data(k).nMice,1)];
    labels(k)          = sprintf('%s (N=%d)', data(k).day, data(k).nMice);
end

%% 4. do boxplot
figure('Color','w', 'Position', [100 100 900 500]);
hBox = boxplot(allVals, allGroup, ...
               'Labels', labels, ...
               'Widths', 0.6, 'Whisker', 1.5);
set(hBox, {'LineWidth'}, {2});

% colour
colors = {[1.0 0.71 0.76],[0.47 0.87 0.47]};
hPatch = findobj(gca, 'Tag', 'Box');
if numel(hPatch) > numel(colors)
    error('colour number less than boxplot number！');
end
for k = 1:numel(hPatch)
    patch(get(hPatch(k), 'XData'), get(hPatch(k), 'YData'), ...
          colors{k}, 'FaceAlpha', 0.5, 'EdgeColor', colors{k});
end

%{
% scatter
hold on
for k = 1:numel(data)
    x = k + (rand(data(k).nMice,1)-0.5)*0.25;
    scatter(x, data(k).responsive, 40, 'MarkerFaceColor',[0.3 0.3 0.3], ...
            'MarkerEdgeColor','k', 'LineWidth',0.5);
end

% mean dashed line
means = cellfun(@mean, {data.responsive});
plot(1:numel(means), means, 'r--', 'LineWidth',1.5);
%}

% axies and title
ylabel('Responsiveness percentage (%)');
xlabel('Groups');
title('Responsive RGC Percentage Across Days');
ylim([0 ceil(max(allVals))]); grid on
set(gca, 'FontSize',12, 'LineWidth',1.2);

% save
saveas(gcf, fullfile(baseFolder, 'Responsive_percentage_boxplot.png'));
close(gcf);

%% 5. stats output
fprintf('\nResponsive Percentage Statistics:\n');
fprintf('%-8s %-4s %-6s %-6s %-6s\n', 'Day','n','Mean','Median','Std');
for k = 1:numel(data)
    resp = data(k).responsive;
    fprintf('%-8d %-4d %-6.1f %-6.1f %-6.1f\n', ...
        data(k).day, numel(resp), mean(resp), median(resp), std(resp));
end
end

