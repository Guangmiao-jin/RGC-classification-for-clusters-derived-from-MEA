function [neuronsprocessed] = processRetinaFlashStimDataFN1_ground_truth(clusterFilepath, stimTime)
spikeTotalLower = 150;
%% get the appropriate paths
filepathPrefix = extractBefore(clusterFilepath, '_cluster');
stimOnFile = [filepathPrefix '_triggerON.npy'];
stimOffFile = [filepathPrefix '_triggerOFF.npy'];
%% load data
matSavePath = extractBefore(clusterFilepath, '.');
matSavePath = [matSavePath '.mat'];

if ~exist(matSavePath)
    data = readHS2_FLAME(clusterFilepath);
    save(matSavePath,"data","-v7.3","-nocompression");
else
    load(matSavePath);
end

stimOnFrames = double(readNPY(stimOnFile));
stimOffFrames = double(readNPY(stimOffFile));
%% split on frames into blocks
blockLimit = 10 * data.Sampling; % 10s


diffOn = diff(stimOnFrames);

% first stim on
blockOnStarts = 1;

% block stim on starts
stimOnBreaks = [blockOnStarts; find(diffOn > blockLimit)+1];

% block stim on ends
stimOnStopBreaks = [ find(diffOn > blockLimit) ;length(diffOn)+1];


for i =1:length(stimOnBreaks)
    stimOnPerBlock{i,:} = stimOnFrames(stimOnBreaks(i):stimOnStopBreaks(i));
end

%% split off frames into blocks

diffOff = diff(stimOffFrames);

% first stim on
blockOffStarts = 1;

% block stim on starts
stimOffBreaks = [blockOffStarts; find(diffOff > blockLimit)+1];

% block stim on ends
stimOffStopBreaks = [ find(diffOff > blockLimit) ;length(diffOff)+1];


for i =1:length(stimOffBreaks)
    stimOffPerBlock{i,:} = stimOffFrames(stimOffBreaks(i):stimOffStopBreaks(i));
end

if length(stimOnPerBlock) ~= 1
    stimOnPerBlock = stimOnPerBlock(1:3);
end
%% start the plotting

cellRasterFolder = extractBefore(clusterFilepath, '.');

if ~exist([cellRasterFolder '_neurons_waveforms'])
    mkdir([cellRasterFolder '_neurons_waveforms']);
end


binsize = 20; % ms
numValidClusters = sum([data.channelNames{6,:}]>spikeTotalLower);

set(0,'DefaultFigureVisible','on');
count = 0;
fs = data.Sampling;

for i = 1:length(data.spiketimestamps)
        if data.channelNames{6,i} >= spikeTotalLower
            spikeTimes = data.spiketimestamps{i};
            if length(stimOnPerBlock) ~= 1
                [~, ~, BI{i}, ISI_on{i}, ISI_off{i}, Ti_on{i}, Ti_off{i}] = find_On_Off_Response(spikeTimes, binsize, fs, stimOnPerBlock(1:3), stimOffPerBlock(1:3), 0.5, stimTime);
            else
                [~, ~, BI{i}, ISI_on{i}, ISI_off{i}, Ti_on{i}, Ti_off{i}] = find_On_Off_Response(spikeTimes, binsize, fs, stimOnPerBlock, stimOffPerBlock, 0.5, stimTime);  
            end
        end
end


num_notRespond = [];
num_OnNeurons.trans = [];
num_OnNeurons.trans_sus = [];
num_OnNeurons.sus = [];
num_OffNeurons.trans = [];
num_OffNeurons.trans_sus = [];
num_OffNeurons.sus = [];
num_OnOffNeurons = [];

for i = 1:length(data.spiketimestamps)
    if data.channelNames{6,i} >= spikeTotalLower
            count = count +1;
            disp(['Analyzing ' num2str(count) ' of ' num2str(numValidClusters)]);

            waveforms = data.waveformsPerCluster{i};
            meanwaveforms = data.waveformClusterMeans(i,:);
            [invertedPeaks,~] = findpeaks(-meanwaveforms);
            numNegativeVally = sum(invertedPeaks>0); % find the neurons containing only one vally of absolute amplitude 1000uV and one peak

            spikenums = data.channelNames{6,i};
            saveName = sprintf('%s/cluster%04d.jpg', [cellRasterFolder '_neurons_waveforms'], i-1);
            plot_waveform_isi(waveforms, meanwaveforms, ISI_on{i}, ISI_off{i}, i, spikenums, BI{i}, Ti_on{i}, Ti_off{i}, saveName);
            
        if numNegativeVally == 1 && abs(invertedPeaks(1)) > 1000
           f = figure('Name', ['Classify Cluster ' num2str(i)], 'Position', [100, 100, 2000, 900], ...
               'Visible', 'on', 'NumberTitle', 'off', 'Color', 'w');

           tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

           nexttile(1);
           img1 = imread(saveName);
           imshow(img1);
           title('Waveform + ISI', 'FontSize', 14, 'FontWeight', 'bold');

           nexttile(2);
           psthName = sprintf('%s/cluster%04d.png', [cellRasterFolder '_PSTHPlots'], i-1);
           if isfile(psthName)
               img2 = imread(psthName);
               imshow(img2);
               title('PSTH Plot', 'FontSize', 14, 'FontWeight', 'bold');
           else
               text(0.3, 0.5, 'PSTH Plot not found', 'FontSize', 16);
               axis off;
            end

           drawnow;
           % === Menu for Classification ==
           choice = menu(['Classify cluster ' num2str(i)], 'ON-trans', 'ON-sus', 'ON-mix', 'OFF-trans', 'OFF-sus', 'OFF-mix', ...
                          'ON-OFF', 'Not responding neuron', 'Not a neuron');           
           close(f);  % Close after classification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            user_labels = {'ON-trans','ON-sus','ON-mix', ...
                           'OFF-trans','OFF-sus','OFF-mix', ...
                           'ON-OFF','Not responding neuron','Not a neuron'};
            user_choice = user_labels{choice};

            % Group assignment
            switch user_choice
                case 'ON-trans'
                    num_OnNeurons.trans = [num_OnNeurons.trans; data.channelNames{4,i}];
                case 'ON-sus'
                    num_OnNeurons.sus = [num_OnNeurons.sus; data.channelNames{4,i}];
                case 'ON-mix'
                    num_OnNeurons.trans_sus = [num_OnNeurons.trans_sus; data.channelNames{4,i}];
                case 'OFF-trans'
                    num_OffNeurons.trans = [num_OffNeurons.trans; data.channelNames{4,i}]; 
                case 'OFF-sus'
                    num_OffNeurons.sus = [num_OffNeurons.sus; data.channelNames{4,i}];
                case 'OFF-mix'
                    num_OffNeurons.trans_sus = [num_OffNeurons.trans_sus; data.channelNames{4,i}];
                case 'ON-OFF'
                    num_OnOffNeurons = [num_OnOffNeurons; data.channelNames{4,i}];
                case 'Not responding neuron'
                    num_notRespond = [num_notRespond; data.channelNames{4,i}];
                case 'Not a neuron'
                    disp('Not a neuron');
            end
            close all;
         else
             label = 'Not a neuron';
        end
    end
end

% Final summary
neuronsprocessed.num_allNeurons = sort([num_OnNeurons.trans;num_OnNeurons.trans_sus; num_OnNeurons.sus; num_OffNeurons.trans; num_OffNeurons.trans_sus; num_OffNeurons.sus;num_OnOffNeurons; num_notRespond]);
neuronsprocessed.num_OnNeurons = num_OnNeurons;
neuronsprocessed.num_OffNeurons = num_OffNeurons;
neuronsprocessed.num_OnOffNeurons = num_OnOffNeurons;
neuronsprocessed.num_notRespond = num_notRespond;


end
