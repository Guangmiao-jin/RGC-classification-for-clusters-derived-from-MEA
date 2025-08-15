function plotNeuronClassificationPie(matfile)
    % Load data
    S = load(matfile);
    S1 = S.totalneurons;
    % Extract directory for saving
    [filepath, name, ~] = fileparts(matfile);

    % Category calculation
    ON_count = 0;
    if isfield(S1.num_OnNeurons, 'trans')
        ON_count = ON_count + length(S1.num_OnNeurons.trans);
    end
    if isfield(S1.num_OnNeurons, 'sus')
        ON_count = ON_count + length(S1.num_OnNeurons.sus);
    end

    OFF_count = 0;
    if isfield(S1.num_OffNeurons, 'trans')
        OFF_count = OFF_count + length(S1.num_OffNeurons.trans);
    end
    if isfield(S1.num_OffNeurons, 'sus')
        OFF_count = OFF_count + length(S1.num_OffNeurons.sus);
    end

    ONOFF_count = 0;
    if isfield(S1, 'num_OnOffNeurons')
        ONOFF_count = length(S1.num_OnOffNeurons);
    end

    unconventional_count = 0;
    if isfield(S1, 'num_unconventional')
        unconventional_count = length(S1.num_unconventional);
    end
    
    if isfield(S1, 'num_notRespond')
        unresponsive_count = length(S1.num_notRespond);
    else
        unresponsive_count = 0;
    end

    % Data for pie chart
    counts = [ON_count, OFF_count, ONOFF_count, unconventional_count, unresponsive_count];
    labels = {'ON', 'OFF', 'ON-OFF', 'Unconventional','Unresponsive'};

    % Remove zero counts for better visualization
    nonzero_idx = counts > 0;
    counts = counts(nonzero_idx);
    labels = labels(nonzero_idx);

    pct       = 100 * counts / sum(counts);
    pieLabels = arrayfun(@(L,P) sprintf('%s\n%.1f%%', L{1}, P), labels, pct, 'UniformOutput', false);


    % Find percentage of transient, sustained and intermediate types in
    % ON-RGCs and OFF-RGCs
    onTrans    = isfield(S1.num_OnNeurons,'trans')    * length(S1.num_OnNeurons.trans);
    onSus      = isfield(S1.num_OnNeurons,'sus')      * length(S1.num_OnNeurons.sus);
    onCounts   = [onTrans, onSus];
    onLabels   = {'Transient','Sustained'};
    onPerc     = 100 * onCounts / sum(onCounts);
    onPieLbls  = arrayfun(@(i) sprintf('%s\n%.1f%%', onLabels{i}, onPerc(i)), 1:length(onLabels), 'UniformOutput',false);


    offTrans   = isfield(S1.num_OffNeurons,'trans')    * length(S1.num_OffNeurons.trans);
    offSus     = isfield(S1.num_OffNeurons,'sus')      * length(S1.num_OffNeurons.sus);
    offCounts  = [offTrans, offSus];
    offLabels  = {'Transient','Sustained'};
    offPerc    = 100 * offCounts / sum(offCounts);
    offPieLbls = arrayfun(@(i) sprintf('%s\n%.1f%%', offLabels{i}, offPerc(i)), 1:length(offLabels), 'UniformOutput',false);

    if sum(onCounts)>0 && sum(offCounts)>0
    figure;
    % ON-RGC pie
    subplot(1,2,1)
    hOn = stablePie(onCounts, onPieLbls);
    title('ON-RGCs');
    axis equal

    % OFF-RGC pie
    subplot(1,2,2)
    hOff = stablePie(offCounts, offPieLbls);
    title('OFF-RGCs');
    axis equal

    % Save
    saveas(gcf, fullfile(filepath, [name '_ON_OFF_Breakdown.png']));
    close(gcf);

    end

    % Create figure with better label handling
    figure;
    h = pie(pct, pieLabels);
    txt = findobj(h, 'Type', 'text');
    for i = 1:numel(txt)
    pos = get(txt(i), 'Position');
    angle = atan2(pos(2), pos(1));
    if i == 4
        r =1.4;
    elseif i == 3
        r =1;
    else
        r =0.7;
    end
    set(txt(i), 'Position', r*[cos(angle), sin(angle), 0], ...  % 
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
    end
    
    for i = 1:numel(txt)
    pos = get(txt(i), 'Position');
    line([0.9*pos(1), pos(1)], [0.9*pos(2), pos(2)], ...
         'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    end

    
    % Set colors
    colorMap = [
        1, 0, 0;        % Red for ON
        0, 0, 1;        % Blue for OFF
        0, 1, 0;        % Green for unconventional
        1, 0.5, 0;      % Orange for ON-OFF
        0.5, 0.5, 0.5   % Gray for Unresponsive
    ];
    colorMap = colorMap(nonzero_idx, :);

    % Apply colors to patches
    patchIdx = 1;
    for k = 1:2:length(h)
        if isfield(h(k), 'FaceColor')
            h(k).FaceColor = colorMap(patchIdx, :);
        elseif isprop(h(k), 'FaceColor')
            set(h(k), 'FaceColor', colorMap(patchIdx, :));
        end
        patchIdx = patchIdx + 1;
    end

    % Save
    saveas(gcf, fullfile(filepath, [name '_NeuronClassificationPie.png']));
    close(gcf);
end

function C = safeField(S, fld)
    if isfield(S, fld)
        C = S.(fld);
    else
        C = {};
    end
end

function hPie = stablePie(counts, labels)
    % get default colour
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    
    % filter out nonzero index
    validIdx = counts > 0;
    counts = counts(validIdx);
    labels = labels(validIdx);
    
    % do the pie plotting
    hPie = pie(counts, labels);
    
    % apply colours in a specific order
    colorIndices = find(validIdx); % get the non-zero position
    for k = 1:numel(counts)
        colorIdx = mod(colorIndices(k)-1, size(defaultColors,1)) + 1;
        set(hPie(2*k-1), 'FaceColor', defaultColors(colorIdx,:));
    end
end

