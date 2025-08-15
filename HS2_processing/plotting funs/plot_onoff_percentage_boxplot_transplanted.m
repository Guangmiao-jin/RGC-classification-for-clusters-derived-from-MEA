function plot_onoff_percentage_boxplot_transplanted(baseFolder)

%% 
groupFolders = dir(fullfile(baseFolder)); 
groupFolders = groupFolders([groupFolders.isdir] & ~ismember({groupFolders.name},{'.','..'}));
data       = [];     % 
id         = 0;

for i = 1:numel(groupFolders)
    dName = groupFolders(i).name;
    dPath = fullfile(baseFolder,dName);
    names = groupFolders(i).name;

    mats = dir(fullfile(dPath,'*_totalneuronsV2.mat'));
    if numel(mats) < 3,  continue,  end   

    onPerc  = nan(numel(mats),1);
    offPerc = nan(numel(mats),1);
    unconPerc = nan(numel(mats),1);

    for j = 1:numel(mats)
        M = load(fullfile(dPath,mats(j).name));
        if ~isfield(M,'totalneurons')
            error('%s no totalneurons struct',mats(j).name);
        end
        T  = M.totalneurons;

        if ~isfield(T,'allneurons')
            error('%s no allneurons struct',mats(j).name);
        end
        totalN = numel(T.allneurons);

        nOn  = 0;
        if isfield(T,'num_OnNeurons')
            nOn = numel(T.num_OnNeurons.trans) + numel(T.num_OnNeurons.sus);
        end

        nOff = 0;
        if isfield(T,'num_OffNeurons')
            nOff = numel(T.num_OffNeurons.trans) + numel(T.num_OffNeurons.sus);
        end

        nUncon = 0;
        if isfield(T,'num_unconventional')
            nUncon = numel(T.num_unconventional);
        end

        onPerc(j)  = nOn  / totalN * 100;
        offPerc(j) = nOff / totalN * 100;
        unconPerc(j) = nUncon / totalN * 100;
    end

    id = id + 1;
    data(id).day   = names;
    data(id).on    = onPerc;
    data(id).off   = offPerc;
    data(id).uncon = unconPerc;
    data(id).nMice = numel(onPerc);
end

if isempty(data)
    disp('⚠️  None of day n sample number ≥ 3，not able to plot。');
    return
end


%% 
Y      = [];      % all values
gDay   = [];      % day n（1,2,3…）
gType  = [];      % 1=ON, 2=OFF
labels = strings(numel(data),1);
labels = arrayfun(@(k) sprintf('%s (N=%d)', ...
    strrep(data(k).day, '_', ' '), ...
    data(k).nMice), ...
    1:numel(data), 'UniformOutput', false);


for k = 1:numel(data)
    nOn  = numel(data(k).on);
    nOff = numel(data(k).off);
    nUncon = numel(data(k).uncon);

    Y      = [Y; data(k).on(:);  data(k).off(:); data(k).uncon(:)];
    gDay   = [gDay; repmat(k, nOn+nOff+nUncon, 1)];
    gType  = [gType; ones(nOn,1); 2*ones(nOff,1);3*ones(nUncon,1)];   % 1/2
end

%% 
figure('Color','w','Position',[100 100 1000 500]);

boxplot(Y, {gDay gType}, ...
        'FactorSeparator',1, ...       % group by gDay, subgroup by gType 
        'Widths',0.6,'Whisker',1.5, ...
        'LabelVerbosity','majorminor');   % only show day n labels

% colours of the boxes：even number (OFF) blue，odd number (ON) red
hBox = flipud(findobj(gca,'Tag','Box'));
for b = 1:numel(hBox)
    if mod(b,3)==0            % reverse the order
        patch(get(hBox(b),'XData'), get(hBox(b),'YData'), [1 0.6 0], 'FaceAlpha',0.4); % unconventional=orange
    elseif mod(b,2)==0
        patch(get(hBox(b),'XData'), get(hBox(b),'YData'), [0 0 1], 'FaceAlpha',0.4); % OFF=blue
    else
        patch(get(hBox(b),'XData'), get(hBox(b),'YData'), [1 0 0], 'FaceAlpha',0.4); % ON=red
    end
end

% adjust the xtick：put XTick to the centre of the two boxes
set(gca,'XTick', 1.5:3.5:4*numel(labels), 'XTickLabel', labels);

ylabel('Responsive neurons (%)');
xlabel('Age');
title('ON / OFF / conventional RGC Responsiveness For transplanted and SHAM injected mice');
ylim([0 max(Y)+1]); grid on
set(gca,'FontSize',12,'LineWidth',1.2);
legend({'ON','OFF', 'unconventional'},'Location','bestoutside');

saveas(gcf, fullfile(baseFolder,'ONOFF_percentage_boxplot.png'));
close(gcf);

%% 
fprintf('\nON / OFF Percentage Statistics\n');
fprintf('%-8s %-4s %-9s %-9s | %-9s %-9s\n',...
        'Day','n','ON-mean','ON-med','OFF-mean','OFF-med');
for k = 1:numel(data)
    fprintf('%-8d %-4d %-9.1f %-9.1f | %-9.1f %-9.1f\n', ...
        data(k).day, data(k).nMice, ...
        mean(data(k).on),  median(data(k).on), ...
        mean(data(k).off), median(data(k).off));
end
end
