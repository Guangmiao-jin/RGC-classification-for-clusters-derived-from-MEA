function plot_onoff_percentage_boxplot(baseFolder)

%% 
dayFolders = dir(fullfile(baseFolder,'day*'));
data       = [];     % 
id         = 0;

for i = 1:numel(dayFolders)
    dName = dayFolders(i).name;
    dPath = fullfile(baseFolder,dName);
    dayNo = str2double(regexp(dName,'\d+','match','once'));

    mats = dir(fullfile(dPath,'*_totalneuronsV2.mat'));
    if numel(mats) < 2,  continue,  end   

    onPerc  = nan(numel(mats),1);
    offPerc = nan(numel(mats),1);

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

        onPerc(j)  = nOn  / totalN * 100;
        offPerc(j) = nOff / totalN * 100;
    end

    id = id + 1;
    data(id).day   = dayNo;
    data(id).on    = onPerc;
    data(id).off   = offPerc;
    data(id).nMice = numel(onPerc);
end

if isempty(data)
    disp('⚠️  None of day n sample number ≥ 3，not able to plot。');
    return
end

% ascending order: day n -> day n+1
[~,ord] = sort([data.day]);
data    = data(ord);

%% 
Y      = [];      % all values
gDay   = [];      % day n（1,2,3…）
gType  = [];      % 1=ON, 2=OFF
labels = strings(numel(data),1);

for k = 1:numel(data)
    nOn  = numel(data(k).on);
    nOff = numel(data(k).off);

    Y      = [Y; data(k).on(:);  data(k).off(:)];
    gDay   = [gDay; repmat(k, nOn+nOff, 1)];
    gType  = [gType; ones(nOn,1); 2*ones(nOff,1)];   % 1/2
    labels(k) = sprintf('Day %d (N=%d)', data(k).day, data(k).nMice);
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
    if mod(b,2)==1            % reverse the order
        patch(get(hBox(b),'XData'), get(hBox(b),'YData'), [1 0 0], 'FaceAlpha',0.4); % ON=red
    else
        patch(get(hBox(b),'XData'), get(hBox(b),'YData'), [0 0 1], 'FaceAlpha',0.4); % OFF=blue
    end
end

% adjust the xtick：put XTick to the centre of the two boxes
set(gca,'XTick', 1.5:2:2*numel(labels), 'XTickLabel', labels);

ylabel('Responsive neurons (%)');
xlabel('Age');
title('ON / OFF RGC Responsiveness Across Days');
ylim([0 max(Y)+5]); grid on
set(gca,'FontSize',12,'LineWidth',1.2);
legend({'ON','OFF'},'Location','bestoutside');

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
