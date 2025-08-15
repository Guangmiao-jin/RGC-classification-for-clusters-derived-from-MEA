function plot_totalneurons_stacked_all(baseFolder)
% set file directory
dayFolders = dir(fullfile(baseFolder, 'day*')); % get all files starting with day

% initialize the data structure
data = struct();
mainCategories = {'OnNeurons', 'OffNeurons', 'OnOffNeurons', 'unconventional', 'notRespond'};
subCategories = {'trans', 'sus', 'trans_sus'};

% iterate through all the folders
for i = 1:length(dayFolders)
    dayName = dayFolders(i).name;
    dayPath = fullfile(baseFolder, dayName);
    
    % extract the number after 'day'
    dayNum = str2double(regexp(dayName, '\d+', 'match'));
    
    % initialize the data structure for storage
    mainData = zeros(1, length(mainCategories));
    onSubData = zeros(1, length(subCategories));
    offSubData = zeros(1, length(subCategories));
    
    % extract all the .mat files in each folder
    matFiles = dir(fullfile(dayPath, '*_totalneuronsV2.mat'));
    
    % iterate through all the .mat files
    for j = 1:length(matFiles)
        matPath = fullfile(dayPath, matFiles(j).name);
        matData = load(matPath);
        matData1 = matData.totalneurons;

        % deal with specific subgroups
        for k = 1:length(mainCategories)
            fieldName = ['num_' mainCategories{k}];
            if isfield(matData1, fieldName)
                if isstruct(matData1.(fieldName)) && ismember(mainCategories{k}, {'OnNeurons', 'OffNeurons'})
                    % if contains on/off subtypes
                    mainData(k) = mainData(k) + length(struct2array(matData1.(fieldName)));
                else
                    % other types
                    mainData(k) = mainData(k) + length(matData1.(fieldName));
                end
            end
        end
        
        % deal with the data fot on and off neurons
        if isfield(matData1, 'num_OnNeurons') && isstruct(matData1.num_OnNeurons)
            for k = 1:length(subCategories)
                if isfield(matData1.num_OnNeurons, subCategories{k})
                    onSubData(k) = onSubData(k) + length(matData1.num_OnNeurons.(subCategories{k}));
                end
            end
        end
        
        if isfield(matData1, 'num_OffNeurons') && isstruct(matData1.num_OffNeurons)
            for k = 1:length(subCategories)
                if isfield(matData1.num_OffNeurons, subCategories{k})
                    offSubData(k) = offSubData(k) + length(matData1.num_OffNeurons.(subCategories{k}));
                end
            end
        end
    end
    
    % store data for each day
    data(i).day = dayNum;
    data(i).mainValues = mainData;
    data(i).onSubValues = onSubData;
    data(i).offSubValues = offSubData;
    x(i) = 10*i;
end

% sort by days
[~, sortIdx] = sort([data.day]);
data = data(sortIdx);

% prepare for plotting!
days = [data.day];
mainY1 = vertcat(data.mainValues)';
a = mainY1;
assignin('base','neurontypes', a);
for i = 1:size(mainY1,2)
    for ii = 1:size(mainY1,1)
        mainY(ii,i) = (mainY1(ii,i) / sum(mainY1(:,i)))*100;
    end
end
b = mainY;
assignin('base','neuronpercs', b);
onSubY1 = vertcat(data.onSubValues)';
offSubY1 = vertcat(data.offSubValues)';
for i = 1:size(onSubY1,2)
    for ii = 1:size(onSubY1,1)
        onSubY(ii,i) = (onSubY1(ii,i) / sum(onSubY1(:,i)))*100;
        offSubY(ii,i) = (offSubY1(ii,i) / sum(offSubY1(:,i)))*100;
    end
end


%% plot first graph and label them
figure(1);
h1 = bar(x, mainY', 'stacked', 'BarWidth', 0.7);

title('Neuron Type Distribution Across Days');
xlabel('Age');
ylabel('Percentage of Neurons (%)');
legend(mainCategories, 'Location', 'bestoutside', 'Interpreter', 'none');
grid on;

% colours
colors1 = lines(length(mainCategories));
for i = 1:length(h1)
    h1(i).FaceColor = colors1(i,:);
end

xlim([min(x)-10, max(x)+10]);
xticks(x);
xticklabels(arrayfun(@(d) sprintf('Day %d', d), days, 'UniformOutput', false));
ylim([0, 100])

% add labels

for i = 1:length(days)
    totalHeight = sum(mainY(:,i));
    cumHeight = 0;
    for j = 1:size(mainY,1)
            cumHeight = cumHeight + mainY(j,i)/2;
            text(x(i), cumHeight, sprintf('%.1f%%', mainY(j,i)),...
                'Color', 'black', 'FontWeight', 'bold',...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
            cumHeight = cumHeight + mainY(j,i)/2;
    end
end
set(gcf, 'Position', [100 100 800 600]);

saveas(gcf, [baseFolder '\Neuron_types_across_days.png'])
close(gcf);

%% second plot
figure(2);

subplot(2,1,1);
ylim([0,100]); hold on;
h2_on = bar(x, onSubY', 'stacked', 'BarWidth', 0.7);
title('On Neuron Subtypes (Transient/Sustained)');
ylabel('Number of Neurons');
legend({'transient','sustained'}, 'Location', 'bestoutside', 'Interpreter', 'none');
grid on;

% Off neurons
subplot(2,1,2);
ylim([0,100]); hold on;
h2_off = bar(x, offSubY', 'stacked', 'BarWidth', 0.7);
title('Off Neuron Subtypes (Transient/Sustained)');
xlabel('Day');
ylabel('Number of Neurons');
legend({'transient','sustained'}, 'Location', 'bestoutside', 'Interpreter', 'none');
grid on; hold on;

% color setting
colors2 = parula(length(subCategories));
for i = 1:length(h2_on)
    h2_on(i).FaceColor = colors2(i,:);
    h2_off(i).FaceColor = colors2(i,:);
end

for i = 1:2
    subplot(2,1,i);
    xlim([min(x)-10, max(x)+10]);
    xticks(x);
    xticklabels(arrayfun(@(d) sprintf('Day %d', d), days, 'UniformOutput', false));
    
    currentAxes = gca;
    bars = findobj(currentAxes, 'Type', 'bar');
    for j = 1:length(days)
        if i == 1
            barData = onSubY(:,j);
        else
            barData = offSubY(:,j);
        end
       
        cumHeight = 0;
        for k = 1:length(barData)
            if barData(k) > 0
                cumHeight = cumHeight + barData(k)/2;
                text(x(j), cumHeight,  sprintf('%.1f%%', barData(k)),...
                    'Color', 'black', 'FontWeight', 'bold',...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
                cumHeight = cumHeight + barData(k)/2;
            end
        end
    end
end


set(gcf, 'Position', [100 100 800 600]);

saveas(gcf, [baseFolder '\Transiency_across_days.png'])
close(gcf);


end