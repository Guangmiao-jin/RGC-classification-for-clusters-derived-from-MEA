function processRetinaFlashStimDataFN1_select(totalneuronsFile)

% Load saved cluster IDs and existing classification
load(totalneuronsFile, 'totalneurons');
baseName = extractBefore(totalneuronsFile,'_totalneurons.mat');
waveformFolder = [baseName '_neurons_waveforms'];
psthFolder     = [baseName '_PSTHPlots'];
progressFile   = [baseName '_classified_temp.mat'];

% Initialize or load existing classifications
if isfile(progressFile)
    load(progressFile);
    fprintf('Loaded previous classification progress.\n');
else
    num_OnNeurons.trans = totalneurons.num_OnNeurons.trans;
    num_OnNeurons.trans_sus = totalneurons.num_OnNeurons.trans_sus;
    num_OnNeurons.sus = totalneurons.num_OnNeurons.sus;
    num_OffNeurons.trans = totalneurons.num_OffNeurons.trans;
    num_OffNeurons.trans_sus = totalneurons.num_OffNeurons.trans_sus;
    num_OffNeurons.sus = totalneurons.num_OffNeurons.sus;
    num_OnOffNeurons.trans_sus = totalneurons.num_OnOffNeurons.trans_sus;
    num_OnOffNeurons.con = totalneurons.num_OnOffNeurons.con;
    num_notRespond = totalneurons.num_notRespond(:);  % Ensure column
    num_notNeuron = [];
    undoStack = {};
end

% Loop through each neuron
totalneurons.raw_all = totalneurons.tobelabelled;
for idx = 1:length(totalneurons.raw_all)
    cluster_id = totalneurons.raw_all(idx);

    % Skip if already classified
    if ismember(cluster_id, [
        ensureColumnNumeric(num_OnNeurons.trans);
        ensureColumnNumeric(num_OnNeurons.trans_sus);
        ensureColumnNumeric(num_OnNeurons.sus);
        ensureColumnNumeric(num_OffNeurons.trans);
        ensureColumnNumeric(num_OffNeurons.trans_sus);
        ensureColumnNumeric(num_OffNeurons.sus);
        ensureColumnNumeric(num_OnOffNeurons.con);
        ensureColumnNumeric(num_OnOffNeurons.trans_sus);
        ensureColumnNumeric(num_notRespond);
        ensureColumnNumeric(num_notNeuron)])
        continue;
    end

    % Load images
    imgName = sprintf('cluster%04d', cluster_id);
    jpgFile = fullfile(waveformFolder, [imgName '.jpg']);
    pngFile = fullfile(psthFolder, [imgName '.png']);

    % Display GUI
    fig = figure(99); clf;
    set(fig, 'Name', sprintf('Neuron %d/%d — Cluster %d', idx, length(totalneurons.raw_all), cluster_id), ...
        'NumberTitle', 'off', 'Position', [100 100 1400 700]);

    try
        subplot(1,2,1); imshow(imread(jpgFile)); title('Waveform');
        subplot(1,2,2); imshow(imread(pngFile)); title('PSTH');
    catch
        warning('Error loading image for Cluster %d. Skipping...', cluster_id);
        continue;
    end

    % Add buttons
    btnLabels = {
        'on_transient', 'on_sus', 'on_trans_sus', ...
        'off_transient', 'off_sus', 'off_trans_sus', ...
        'on_off_con', 'on_off_trans_sus', ...
        'not_responding', 'not_a_neuron'
    };

    assignments = {
        'num_OnNeurons.trans', 'num_OnNeurons.sus', 'num_OnNeurons.trans_sus', ...
        'num_OffNeurons.trans', 'num_OffNeurons.sus', 'num_OffNeurons.trans_sus', ...
        'num_OnOffNeurons.con', 'num_OnOffNeurons.trans_sus', ...
        'num_notRespond', 'num_notNeuron'
    };

    for i = 1:length(btnLabels)
        uicontrol('Style', 'pushbutton', 'String', btnLabels{i}, ...
            'FontSize', 12, 'Position', [50 + (i-1)*120, 20, 110, 40], ...
            'Callback', @(src, event) classify_and_continue(assignments{i}));
    end

    % Add undo button
    uicontrol('Style', 'pushbutton', 'String', 'Undo (u)', ...
        'FontSize', 12, 'BackgroundColor', [1 0.5 0.5], ...
        'Position', [1250, 20, 100, 40], ...
        'Callback', @(src, event) undo_last());

    % Wait for classification
    uiwait(fig);
end

% Combine all classified clusters
neuronsprocessed.num_allNeurons = sort(vertcat(...
    ensureColumnNumeric(num_OnNeurons.trans), ...
    ensureColumnNumeric(num_OnNeurons.trans_sus), ...
    ensureColumnNumeric(num_OnNeurons.sus), ...
    ensureColumnNumeric(num_OffNeurons.trans), ...
    ensureColumnNumeric(num_OffNeurons.trans_sus), ...
    ensureColumnNumeric(num_OffNeurons.sus), ...
    ensureColumnNumeric(num_OnOffNeurons.con), ...
    ensureColumnNumeric(num_OnOffNeurons.trans_sus), ...
    ensureColumnNumeric(num_notRespond)));

% Save final
save([baseName '_neurons_classified.mat'], 'num_OnNeurons', 'num_OffNeurons', ...
    'num_OnOffNeurons', 'num_notRespond', 'num_notNeuron', 'neuronsprocessed');
fprintf('Classification complete.\n');

    %% === NESTED FUNCTIONS ===

    function classify_and_continue(fieldName)
        dotIdx = strfind(fieldName, '.');

        if isempty(dotIdx)
            % Scalar array case
            switch fieldName
                case 'num_notRespond'
                    num_notRespond = [num_notRespond(:); cluster_id];
                case 'num_notNeuron'
                    num_notNeuron = [num_notNeuron(:); cluster_id];
            end
        else
            structVar = fieldName(1:dotIdx-1);
            field = fieldName(dotIdx+1:end);

            switch structVar
                case 'num_OnNeurons'
                    current = num_OnNeurons.(field);
                    current = [current(:); cluster_id];
                    num_OnNeurons.(field) = current;
                case 'num_OffNeurons'
                    current = num_OffNeurons.(field);
                    current = [current(:); cluster_id];
                    num_OffNeurons.(field) = current;
                case 'num_OnOffNeurons'
                    current = num_OnOffNeurons.(field);
                    current = [current(:); cluster_id];
                    num_OnOffNeurons.(field) = current;
            end
        end

        undoStack{end+1} = struct('cluster_id', cluster_id, 'field', fieldName);
        save(progressFile, 'num_OnNeurons', 'num_OffNeurons', ...
            'num_OnOffNeurons', 'num_notRespond', 'num_notNeuron', 'undoStack');
        disp(['Saved cluster ' num2str(cluster_id) ' to ' fieldName]);
        uiresume(fig);
        close(fig);
    end

    function undo_last()
        if isempty(undoStack)
            disp('Nothing to undo!');
            return;
        end
        last = undoStack{end};
        undoStack(end) = [];

        dotIdx = strfind(last.field, '.');
        if isempty(dotIdx)
            switch last.field
                case 'num_notRespond'
                    num_notRespond = setdiff(num_notRespond, last.cluster_id);
                case 'num_notNeuron'
                    num_notNeuron = setdiff(num_notNeuron, last.cluster_id);
            end
        else
            structVar = last.field(1:dotIdx-1);
            field = last.field(dotIdx+1:end);

            switch structVar
                case 'num_OnNeurons'
                    num_OnNeurons.(field) = setdiff(num_OnNeurons.(field), last.cluster_id);
                case 'num_OffNeurons'
                    num_OffNeurons.(field) = setdiff(num_OffNeurons.(field), last.cluster_id);
                case 'num_OnOffNeurons'
                    num_OnOffNeurons.(field) = setdiff(num_OnOffNeurons.(field), last.cluster_id);
            end
        end

        save(progressFile, 'num_OnNeurons', 'num_OffNeurons', ...
            'num_OnOffNeurons', 'num_notRespond', 'num_notNeuron', 'undoStack');
        disp(['Undid classification of Cluster ' num2str(last.cluster_id)]);
        uiresume(fig);
        close(fig);
    end

    function x = ensureColumnNumeric(x)
        if isempty(x)
            x = [];
        elseif ~isnumeric(x)
            x = [];
        else
            x = x(:);  % force column vector
        end
    end

end
