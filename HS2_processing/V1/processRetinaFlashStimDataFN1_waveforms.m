function processRetinaFlashStimDataFN1_waveforms(clusterFilepath, stimTime, spikeTotalLower)

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
 
totalneurons.raw_all = [];


for i = 1:length(data.spiketimestamps)
    if data.channelNames{6,i} >= spikeTotalLower
            count = count +1;

            waveforms = data.waveformsPerCluster{i};
            meanwaveforms = data.waveformClusterMeans(i,:);
            [invertedPeaks,~] = findpeaks(-meanwaveforms);
            numNegativeVally = sum(invertedPeaks>0); % find the neurons containing only one vally of absolute amplitude 1000uV and one peak

            spikenums = data.channelNames{6,i};
            saveName = sprintf('%s/cluster%04d.jpg', [cellRasterFolder '_neurons_waveforms'], i-1);
            plot_waveform_isi(waveforms, meanwaveforms, ISI_on{i}, ISI_off{i}, i, spikenums, BI{i}, Ti_on{i}, Ti_off{i}, saveName);
            disp(['Analyzing ' num2str(count) ' of ' num2str(numValidClusters)]);

        if numNegativeVally == 1 && abs(invertedPeaks(1)) > 1000
           totalneurons.raw_all = [totalneurons.raw_all; data.channelNames{4,i}];
           disp(['Finishing ' num2str(count) ' of ' num2str(numValidClusters)]);
        end
    end
end

% Final summary
save([cellRasterFolder '_totalneurons.mat'],"totalneurons");

end
