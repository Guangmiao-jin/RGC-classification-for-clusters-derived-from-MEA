function f = plot_waveform_isi_test(waveforms, meanwaveforms, ISI_on, ISI_off, clusterID, spikeCount, BI, TI_on, TI_off, ratio, pvalues, savePath)
    
    % Create figure
    f = figure('Visible', 'off', ...
               'Renderer', 'painters', ...
               'Position', [100 100 900 300]);
    tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

    colors_on = [0.3 0.7 0.9; 0.4 0.9 0.5; 0.6 0.5 0.9];
    colors_off = [0.1 0.3 0.8; 0.1 0.6 0.2; 0.4 0.2 0.6];
    %condition_labels = {'scotopic', 'mesopic', 'photopic'};

    % --- Plot 1: Waveforms
    nexttile(1);
    hold on;
    batchSize = min(100, size(waveforms, 1));
    for i = 1:batchSize:size(waveforms, 1)
        idx = i:min(i+batchSize-1, size(waveforms, 1));
        plot(waveforms(idx,:)', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.1);
    end
    plot(meanwaveforms, 'k', 'LineWidth', 1.5);
    xlabel('Time (samples)');
    title('Waveforms');
    hold off;
    
    % --- Plot 2: ISI ON (all conditions together)
    nexttile(2);
    hold on;
    for cond = 1:numel(ISI_on)
        [counts, edges] = histcounts(ISI_on{cond}, 50, 'Normalization', 'probability');
        stairs(edges(1:end-1), counts, 'Color', colors_on(cond,:), 'LineWidth', 1);
    end
    xlabel('ISI ON (s)');
    title(['TI ON: ' strjoin(arrayfun(@num2str, TI_on, 'UniformOutput', false), '/')]);
    hold off;

    % --- Plot 3: ISI OFF (all conditions together)
    nexttile(3);
    hold on;
    for cond = 1:numel(ISI_off)
        [counts, edges] = histcounts(ISI_off{cond}, 50, 'Normalization', 'probability');
        stairs(edges(1:end-1), counts, 'Color', colors_off(cond,:), 'LineWidth', 1);
    end
    xlabel('ISI OFF (s)');
    title(['TI OFF: ' strjoin(arrayfun(@num2str, TI_off, 'UniformOutput', false), '/')]);
    hold off;

    % ---  Title
    ratio_str = sprintf('Ratio: %s', strjoin(arrayfun(@(x) sprintf('%.2f', x), ratio, 'UniformOutput', false), '/'));
    pval_str = sprintf('p-value: %s', strjoin(arrayfun(@(x) sprintf('%.3f', x), pvalues, 'UniformOutput', false), '/'));

    stats_title = sprintf('%s\n%s', ratio_str, pval_str);
    
    bi_str = sprintf('BI: %s', strjoin(arrayfun(@(x) sprintf('%.2f', x), BI, 'UniformOutput', false), '/'));
    main_title = sprintf('Cluster ID: %d (Spks: %d)\n%s\n%s', ...
        clusterID-1, spikeCount, bi_str, stats_title);
    
    sgtitle(main_title, 'FontSize', 10, 'FontWeight', 'bold');
    % --- Save figure
    exportgraphics(f, savePath, 'Resolution', 150);
    close(f);
    delete(f);
    drawnow;
end
