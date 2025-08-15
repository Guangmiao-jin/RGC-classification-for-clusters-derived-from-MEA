function processRetinaFlashStimDataFN1_waveforms_test(clusterFilepath, stimTime, spikeTotalLower, minspikes)
%% disenable figure display
set(0, 'DefaultFigureVisible', 'off');  
filepathPrefix = extractBefore(clusterFilepath, '_cluster');
matSavePath = strrep(clusterFilepath, '.hdf5', '.mat');  

%% load data
if ~exist(matSavePath, 'file')
    fprintf('Loading raw data...\n');
    data = readHS2_FLAME(clusterFilepath);
    save(matSavePath, "data", "-v7.3", "-nocompression");  % large hdf5 file size > 2GB
else
    fprintf('Loading pre-saved data...\n');
    data = load(matSavePath);
    data = data.data; 
end

%% load stimulation time points
stimOnFile = [filepathPrefix '_triggerON.npy'];
stimOffFile = [filepathPrefix '_triggerOFF.npy'];

% extract npy file
if exist('parreadNPY.m', 'file') == 0
    stimOnFrames = double(readNPY(stimOnFile));
    stimOffFrames = double(readNPY(stimOffFile));
else
    stimOnFrames = double(parreadNPY(stimOnFile));
    stimOffFrames = double(parreadNPY(stimOffFile));
end

%% split time
blockLimit = 10 * data.Sampling;  % 10s
[stimOnPerBlock, stimOffPerBlock] = splitStimBlocks(stimOnFrames, stimOffFrames, blockLimit);

%% Select valid neuron ids
validNeurons = find([data.channelNames{6,:}] >= spikeTotalLower);
numValidClusters = length(validNeurons);

% pre-allocate
results = struct(...
    'BI', cell(1, numValidClusters), ...
    'ISI_on', cell(1, numValidClusters), ...
    'ISI_off', cell(1, numValidClusters), ...
    'Ti_on', cell(1, numValidClusters), ...
    'Ti_off', cell(1, numValidClusters),...
    'ratio', cell(1, numValidClusters),...
    'pvalues', cell(1, numValidClusters));


% processing
for idx = 1:numValidClusters
    i = validNeurons(idx);
    fprintf('Processing neuron %d/%d\n', idx, numValidClusters);
    spikeTimes = data.spiketimestamps{i};
    if length(stimOnPerBlock) ~= 1
        [~, ~, results(idx).BI, results(idx).ISI_on, results(idx).ISI_off, results(idx).Ti_on, results(idx).Ti_off, results(idx).ratio,results(idx).pvalues] = find_On_Off_Response_test(spikeTimes, 20, data.Sampling, stimOnPerBlock(1:3), stimOffPerBlock(1:3), 0.5, stimTime);
    else
        [~, ~, results(idx).BI, results(idx).ISI_on, results(idx).ISI_off, results(idx).Ti_on, results(idx).Ti_off, results(idx).ratio,results(idx).pvalues] = find_On_Off_Response_test(spikeTimes, 20, data.Sampling, stimOnPerBlock, stimOffPerBlock, 0.5, stimTime);
    end
end

%% wavform and plotting
cellRasterFolder = extractBefore(clusterFilepath, '_cluster');
waveformDir = [cellRasterFolder '_cluster_neurons_waveforms'];
if ~exist(waveformDir, 'dir')
    mkdir(waveformDir);
end

%totalneurons.raw_all = [];
totalneurons.num_OnNeurons.trans = [];
totalneurons.num_OnNeurons.sus = [];
totalneurons.num_OnNeurons.trans_sus = [];

totalneurons.num_OffNeurons.trans = [];
totalneurons.num_OffNeurons.sus = [];
totalneurons.num_OffNeurons.trans_sus = [];

totalneurons.num_OnOffNeurons.con = [];
totalneurons.num_OnOffNeurons.trans_sus = [];
totalneurons.tobelabelled = [];

%totalneurons.unknown = [];
totalneurons.num_notRespond = [];

for idx = 1:numValidClusters
    i = validNeurons(idx);
    waveforms = data.waveformsPerCluster{i};
    meanwaveforms = data.waveformClusterMeans(i,:);
    
    % detect negative peaks and positive peaks
    [invertedPeaks, t0] = findpeaks(-meanwaveforms);
    numNegativeVally = sum(invertedPeaks>0); 
    [~,t1] = max(meanwaveforms);
    
    % save those neuron ids which fall into given conditions
if numNegativeVally == 1 && abs(invertedPeaks(1)) > 1000 % make sure there is only one negative peak and absolute value is greater than 1000
    peaktimediff = t1-t0(1);
  if peaktimediff > 2
    if all(results(idx).pvalues < 0.05) && data.channelNames{6,i} > minspikes
                if all(results(idx).BI >=0.33)
                    if mean(results(idx).Ti_on) > 0.85
                         totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,i}; 
                    end
                elseif mean(results(idx).BI) > 0.3
                    if mean(results(idx).Ti_on) > 0.85
                         totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,i}; 
                    end
                elseif all(results(idx).BI <=-0.33)
                    if mean(results(idx).Ti_off) > 0.85
                         totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,i};
                    end   
                elseif mean(results(idx).BI) < -0.3
                    if mean(results(idx).Ti_off) > 0.85
                         totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,i};
                    end
                else
                    totalneurons.tobelabelled(end+1) = data.channelNames{4,i};
                end            
    else
        if data.channelNames{6,i} <= 300
            continue;
        else
            if all(results(idx).BI >=0.33)
                    if mean(results(idx).Ti_on) > 0.85
                         totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,i}; 
                    end
            elseif all(results(idx).BI <=-0.33)
                    if mean(results(idx).Ti_off) > 0.85
                         totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,i}; 
                    else
                        totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,i};
                    end
            elseif all(results(idx).pvalues >= 0.05)
                totalneurons.num_notRespond(end+1) = data.channelNames{4,i};
            else
                totalneurons.tobelabelled(end+1) = data.channelNames{4,i};
           end            
        end
    end
  end
end
    
    % plotting and saving in local directory
    
    saveName = sprintf('%s/cluster%04d.jpg', waveformDir, i-1);
    h = plot_waveform_isi_test(waveforms, meanwaveforms, ...
        results(idx).ISI_on, results(idx).ISI_off, i, ...
        data.channelNames{6,i}, results(idx).BI, ...
        results(idx).Ti_on, results(idx).Ti_off, results(idx).ratio,results(idx).pvalues,saveName);
    
    
    % clearing un-wanted variables
    if mod(idx, 50) == 0
        clear waveforms meanwaveforms invertedPeaks;
    end
end

%% saving final results
save([cellRasterFolder '_cluster_totalneurons.mat'], 'totalneurons');  % 
fprintf('Analysis completed!\n');
end

%% function that separates on and off into blocks
function [onBlocks, offBlocks] = splitStimBlocks(onFrames, offFrames, blockLimit)

    onDiff = diff(onFrames);
    onBreaks = [1; find(onDiff > blockLimit) + 1];
    onStops = [find(onDiff > blockLimit); length(onFrames)];
    onBlocks = arrayfun(@(s,e) onFrames(s:e), onBreaks, onStops, 'UniformOutput', false);
    
    offDiff = diff(offFrames);
    offBreaks = [1; find(offDiff > blockLimit) + 1];
    offStops = [find(offDiff > blockLimit); length(offFrames)];
    offBlocks = arrayfun(@(s,e) offFrames(s:e), offBreaks, offStops, 'UniformOutput', false);
end