function plot_locations_test(clusterFilepath, neuronIDpath)
    %% Enable figure display
    load(neuronIDpath);
    set(0, 'DefaultFigureVisible', 'on'); 
    filepathPrefix = extractBefore(clusterFilepath, '_cluster');
    matSavePath = strrep(clusterFilepath, '.hdf5', '.mat'); 

    %% Load data
    if ~exist(matSavePath, 'file')
        fprintf('Loading raw data...\n');
        data = readHS2_FLAME(clusterFilepath);
        save(matSavePath, "data", "-v7.3", "-nocompression");  % large hdf5 file size > 2GB
    else
        fprintf('Loading pre-saved data...\n');
        data = load(matSavePath);
        data = data.data;
    end

    %% Prepare neuron ID lists
    ON_all = sort(cat(1, num_OnNeurons.trans(:), num_OnNeurons.sus(:), num_OnNeurons.trans_sus(:)));
    OFF_all = sort(cat(1, num_OffNeurons.trans(:), num_OffNeurons.sus(:), num_OffNeurons.trans_sus(:)));
    ON_OFF_all = sort(cat(1, num_OnOffNeurons.con(:), num_OnOffNeurons.trans_sus(:)));
    not_responding = num_notRespond;

    %% Convert channel names to string if needed
    channelNames = string(data.channelNames(4, :));

    %% Match indices
    [~, IDX_on] = ismember(string(ON_all), channelNames);
    [~, IDX_off] = ismember(string(OFF_all), channelNames);
    [~, IDX_on_off] = ismember(string(ON_OFF_all), channelNames);
    [~, IDX_not_respond] = ismember(string(not_responding), channelNames);

    %% Filter valid indices (nonzero)
    IDX_on = IDX_on(IDX_on > 0);
    IDX_off = IDX_off(IDX_off > 0);
    IDX_on_off = IDX_on_off(IDX_on_off > 0);
    IDX_not_respond = IDX_not_respond(IDX_not_respond > 0);

    %% Plotting
    figure;
    hold on;
    grid on;
    xlim([0 3780]);
    ylim([0 3780]);
    title('RGC type distribution map');

    scatter(data.centres(1, IDX_on), data.centres(2, IDX_on), 50, 'r', 'filled');
    scatter(data.centres(1, IDX_off), data.centres(2, IDX_off), 50, 'b', 'filled');
    scatter(data.centres(1, IDX_on_off), data.centres(2, IDX_on_off), 50, 'g', 'filled');
    scatter(data.centres(1, IDX_not_respond), data.centres(2, IDX_not_respond), 50, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',0.5);

    legend('ON-type', 'OFF-type', 'ON-OFF-type', 'Neurons not responding', 'Location', 'best');

    %% Save figure
    outputDir = [filepathPrefix '_cluster_neurons_waveforms'];
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    saveName = fullfile(outputDir, 'scatter_graph.jpg');
    saveas(gcf, saveName);
end
