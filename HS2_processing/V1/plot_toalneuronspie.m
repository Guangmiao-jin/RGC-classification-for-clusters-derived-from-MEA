function plot_toalneuronspie(matfile)

%%%%%%% For V1 processing %%%%%%%%%%
    % Load data (_classified.mat)
    S = load(matfile);
    % Extract directory for saving
    [filepath, name, ~] = fileparts(matfile);

    % Category calculation
    ON_count = 0;
    if isfield(S.num_OnNeurons, 'trans')
        ON_count = ON_count + length(S.num_OnNeurons.trans);
    end
    if isfield(S.num_OnNeurons, 'sus')
        ON_count = ON_count + length(S.num_OnNeurons.sus);
    end

    OFF_count = 0;
    if isfield(S.num_OffNeurons, 'trans')
        OFF_count = OFF_count + length(S.num_OffNeurons.trans);
    end
    if isfield(S.num_OffNeurons, 'sus')
        OFF_count = OFF_count + length(S.num_OffNeurons.sus);
    end

    ONOFF_count = 0;
    if isfield(S.num_OnOffNeurons, 'con')
        ONOFF_count = length(S.num_OnOffNeurons.con);
    end

    unconventional_count = 0;
    if isfield(S.num_OnNeurons, 'trans_sus')
        unconventional_count = unconventional_count + length(S.num_OnNeurons.trans_sus);
    end
    if isfield(S.num_OffNeurons, 'trans_sus')
        unconventional_count = unconventional_count + length(S.num_OffNeurons.trans_sus);
    end
    if isfield(S.num_OnOffNeurons, 'trans_sus')
        unconventional_count = unconventional_count + length(S.num_OnOffNeurons.trans_sus);
    end

    if isfield(S, 'num_notRespond')
        unresponsive_count = length(S.num_notRespond);
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
    onTrans    = isfield(S.num_OnNeurons,'trans')    * length(S.num_OnNeurons.trans);
    onSus      = isfield(S.num_OnNeurons,'sus')      * length(S.num_OnNeurons.sus);
    onInter    = isfield(S.num_OnNeurons,'trans_sus')* length(S.num_OnNeurons.trans_sus);
    onCounts   = [onTrans, onSus, onInter];
    onLabels   = {'Transient','Sustained','Intermediate'};
    onPerc     = 100 * onCounts / sum(onCounts);
    onPieLbls  = arrayfun(@(i) sprintf('%s\n%.1f%%', onLabels{i}, onPerc(i)), 1:3, 'UniformOutput',false);


    offTrans   = isfield(S.num_OffNeurons,'trans')    * length(S.num_OffNeurons.trans);
    offSus     = isfield(S.num_OffNeurons,'sus')      * length(S.num_OffNeurons.sus);
    offInter   = isfield(S.num_OffNeurons,'trans_sus')* length(S.num_OffNeurons.trans_sus);
    offCounts  = [offTrans, offSus, offInter];
    offLabels  = {'Transient','Sustained','Intermediate'};
    offPerc    = 100 * offCounts / sum(offCounts);
    offPieLbls = arrayfun(@(i) sprintf('%s\n%.1f%%', offLabels{i}, offPerc(i)), 1:3, 'UniformOutput',false);

    figure;
    % ON-RGC pie
    subplot(1,2,1)
    hOn = pie(onCounts, onPieLbls);
    title('ON-RGCs');
    axis equal

    % OFF-RGC pie
    subplot(1,2,2)
    hOff = pie(offCounts, offPieLbls);
    title('OFF-RGCs');
    axis equal

    % Save
    saveas(gcf, fullfile(filepath, [name 'V1_ON_OFF_Breakdown.png']));
    close(gcf);


    % Create figure with better label handling
    figure;
    h = pie(pct, pieLabels);
    txt = findobj(h, 'Type', 'text');
    for i = 1:numel(txt)
    pos = get(txt(i), 'Position');
    angle = atan2(pos(2), pos(1));
    if i == 4
        r =1.2;
    elseif i == 3
        r =1;
    else
        r =0.8;
    end
    set(txt(i), 'Position', r*[cos(angle), sin(angle), 0], ...  % 固定距离
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
    saveas(gcf, fullfile(filepath, [name 'V1_NeuronClassificationPie.png']));
    close(gcf);
end

function C = safeField(S, fld)
    if isfield(S, fld)
        C = S.(fld);
    else
        C = {};
    end
end