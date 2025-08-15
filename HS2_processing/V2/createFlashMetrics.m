function [responseMetrics] = createFlashMetrics(data, stimON_Events, stimOFF_Events, ~, ~, prestimTimePSTH, postStimTimePSTH)

% preAlignTime = 1; %s
% postAlignTime = 1; %s
binsize = 30; % ms
fs = data.Sampling;
reclen = max(data.times)/fs;

%% Run through all the clusters
% for each cluster
spikes = data.spiketimestamps;
parfor i = 1:length(data.spiketimestamps)
    spikeFrames = spikes{i};
    %% ISI violation value
    ISI{i} = diff(spikeFrames/fs);
    %ISI_mean(i) = mean(ISI{i});
    num_vio = sum(ISI{i} < 0.001);
    ISI_min = min(ISI{i});
    ISI_vio(i) = abs((num_vio * reclen) / (2*length(spikeFrames)^2*(0.001-ISI_min))); % include violation value less than 5%

    %% Comparison between pre- and post-stimulation PSTH for ON- and OFF-flashes
    [ZvaluePerClusterBlkON{i}, ZvaluePerClusterBlkOFF{i}] = stimAlignedrank(spikeFrames, fs, stimON_Events, stimOFF_Events , postStimTimePSTH);

    %% Median absolute divation?
    %% create trial based PSTHs
    warning('off','MATLAB:colon:operandsNotRealScalar'); % stops Warning: Colon operands must be real scalars. This warning will become an error in a future release.
    [trialPSTHs{i}, trialSpikeStruct{i}, binEdges{i}] = createTrialPSTHs(spikeFrames, fs, stimON_Events, stimOFF_Events, prestimTimePSTH, postStimTimePSTH);
    warning('on','MATLAB:colon:operandsNotRealScalar');

    %% get response quality
    QI(i,:) = retinaResponseQuality(trialPSTHs{i});

    %% bias index, ratio and tau times

    [ISIon_mean{i}, ISIoff_mean{i}, BI{i}, ~, ~, Ti_on{i}, Ti_off{i}, ratio{i}, ~] = find_On_Off_Response_test(spikeFrames, fs, stimON_Events, stimOFF_Events);
end

%% Put everything into struct
responseMetrics.ISI_vio = ISI_vio;
responseMetrics.ISIon_mean = ISIon_mean;
responseMetrics.ISIoff_mean = ISIoff_mean;
responseMetrics.ZValueON = ZvaluePerClusterBlkON;
responseMetrics.ZValueOFF = ZvaluePerClusterBlkOFF;
responseMetrics.trialSpikes = trialSpikeStruct;
responseMetrics.trialPSTHs = trialPSTHs;
responseMetrics.PSTH_binEdges = binEdges{1};
responseMetrics.responseQuality = QI;
responseMetrics.BI = BI;
responseMetrics.TI_on = Ti_on;
responseMetrics.TI_off = Ti_off;
responseMetrics.ratio = ratio;
end