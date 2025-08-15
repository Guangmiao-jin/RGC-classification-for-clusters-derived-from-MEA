function plot_waveform_isi(waveforms, meanwaveforms, ISI_on, ISI_off, clusterID, spikeCount, BI, TI_on, TI_off, savePath)
    numConditions = numel(ISI_on);
    
    % Create figure
    f = figure('Visible', 'off');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Define colors
    colors_on = [0.3 0.7 0.9; 0.4 0.9 0.5; 0.6 0.5 0.9];  % Light colors for ON
    colors_off = [0.1 0.3 0.8; 0.1 0.6 0.2; 0.4 0.2 0.6]; % Dark colors for OFF
    condition_labels = {'scotopic', 'mesopic', 'photopic'};  % Human-readable labels

    % --- Plot 1: Waveforms
    nexttile(1);
    plot(waveforms', 'Color', [0.6 0.6 0.6]); hold on;
    plot(meanwaveforms, 'k', 'LineWidth', 2);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    title('Waveforms', 'FontWeight', 'bold');
    set(gca, 'FontSize', 5);
    hold off;

    % --- Plot 2: ISI ON (all conditions together)
    nexttile(2);
    hold on;
    legend_entries_on = {};
    for cond = 1:numConditions
        histogram(ISI_on{cond}, 'Normalization', 'probability', ...
            'FaceColor', colors_on(cond,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        legend_entries_on{end+1} = [condition_labels{cond} ' ON'];
    end
    hold off;
    xlabel('ISI ON (s)');
    ylabel('Probability');
    ti_on_string_parts = arrayfun(@(x) sprintf('%.2f', x), TI_on, 'UniformOutput', false);
    TI_on_string = ['TI ON: ' strjoin(ti_on_string_parts, '/')];
    title(TI_on_string, 'FontSize', 6);
    legend(legend_entries_on, 'Location', 'best');
    set(gca, 'FontSize', 5);

    % --- Plot 3: ISI OFF (all conditions together)
    nexttile(3);
    hold on;
    legend_entries_off = {};
    for cond = 1:numConditions
        histogram(ISI_off{cond}, 'Normalization', 'probability', ...
            'FaceColor', colors_off(cond,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        legend_entries_off{end+1} = [condition_labels{cond} ' OFF'];
    end
    hold off;
    xlabel('ISI OFF (s)');
    ylabel('Probability');
    ti_off_string_parts = arrayfun(@(x) sprintf('%.2f', x), TI_off, 'UniformOutput', false);
    TI_off_string = ['TI OFF: ' strjoin(ti_off_string_parts, '/')];
    title(TI_off_string, 'FontSize', 6);
    legend(legend_entries_off, 'Location', 'best');
    set(gca, 'FontSize', 5);

    % ---  Title
    bi_string_parts = arrayfun(@(x) sprintf('%.2f', x), BI, 'UniformOutput', false);
    BI_string = ['BI: ' strjoin(bi_string_parts, '/')];
    
    sgtitle(sprintf('Cluster ID: %d (Spks: %d) %s', clusterID-1, spikeCount, BI_string), 'FontSize', 10, 'FontWeight', 'bold');

    % --- Save figure
    exportgraphics(f, savePath, 'Resolution', 150);
    close(f);
end
