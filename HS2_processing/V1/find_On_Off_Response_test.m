function [ISIon_mean, ISIoff_mean, BI, ISI_on, ISI_off, Ti_on, Ti_off, ratio, pvalues] = find_On_Off_Response_test(spikeFrames, fs, stimON_Events, stimOFF_Events)
%% ON stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 30;
postAlignTime = 2;
if isscalar(stimON_Events)
    nBlocks = 1;
else
    nBlocks = 3;
end

for stimblk = 1:nBlocks
    % alignment times in sec
    postAlignInFrames = postAlignTime * fs;
    trSpikesAll= [];
    trSpikesStructON = [];
    % for each trial
    for tr = 1:length(stimON_Events{stimblk})

        % trial times in frames
        trStart = stimON_Events{stimblk}(tr);
        trEnd = stimON_Events{stimblk}(tr)+ postAlignInFrames;
        trialLenSec = (trEnd-trStart)/fs;

        % find inclusions for trial start and end
        trStartIndx = find(spikeFrames >trStart,1, 'first');
        trEndIndx = find(spikeFrames <trEnd,1, 'last');

        % get the spikes
        trSpikes = spikeFrames(trStartIndx:  trEndIndx);
        trSpikes = trSpikes / fs; % convert to sec
        trSpikes = trSpikes- (stimON_Events{stimblk}(tr)/ fs); % rezero to alignment event

        % put into cell array for psth
        trSpikesStructON{tr} = trSpikes;
        trSpikesAll = [trSpikesAll; length(trSpikes)];
    end

    trSpikesMeanOn(stimblk) = sum(trSpikesAll)/length(stimON_Events{stimblk});

    %% PSTH
    trialLen            = mean(trialLenSec) * 1000;                % trial length ms
    edges               = 0: (binsize/1000): trialLenSec;
    binActual           = edges;
    nbins               = trialLen/binsize;                        % Bin duration in [ms]
    nobins              = 1000/binsize;                            % sec / bin
    all = [];
    ISI_on{stimblk} = [];
    trSpikesVarOnAll = [];
    for iTrial = 1:length(trSpikesStructON)
        all             = [all; trSpikesStructON{iTrial}];               % Concatenate spikes of all trials
        ISI_on{stimblk} = [ISI_on{stimblk}; diff(trSpikesStructON{iTrial})];
        trSpikesVarOnAll = [trSpikesVarOnAll; (trSpikesAll(iTrial) - trSpikesMeanOn(stimblk))^2];
        [trialHistTemp, ~] = histcounts(trSpikesStructON{iTrial},binActual);
        trialHistson{stimblk}(iTrial,:) = trialHistTemp;
    end

   
   ISIon_mean(stimblk) = mean(ISI_on{stimblk});
   meanPSTH = mean(trialHistson{stimblk});
   ONcountAverageSec     = (meanPSTH) * nobins;
   [A1(stimblk), peak_bin(stimblk)] = max(ONcountAverageSec);
    A2(stimblk) = max(ONcountAverageSec)*(1/exp(1));
    % Find the index after the peak where the PSTH drops below A2
    post_peak = ONcountAverageSec(peak_bin(stimblk):end);
    below_A2 = find(post_peak < A2(stimblk), 1, 'first');
    if ~isempty(below_A2)
        tau_bin = peak_bin(stimblk) + below_A2 - 1; % Index in original array
        Ti_on(stimblk) = (tau_bin - peak_bin(stimblk)) * (binsize/1000); % Time difference (ms if binsize in ms)
    else
        below_A2 = ONcountAverageSec(end);
        tau_bin = peak_bin(stimblk) + below_A2 - 1; % Index in original array
        Ti_on(stimblk) = (tau_bin - peak_bin(stimblk)) * (binsize/1000); % Time difference (ms if binsize in ms)
    end

end

%% OFF stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stimblk = 1:nBlocks
    % alignment times in sec
    postAlignInFrames = postAlignTime * fs;
    trSpikesAll= [];
    trSpikesStructOFF = [];

    % for each trial
    for tr = 1:length(stimOFF_Events{stimblk})

        % trial times in frames
        trStart = stimOFF_Events{stimblk}(tr);
        trEnd = stimOFF_Events{stimblk}(tr)+ postAlignInFrames;
        trialLenSec = (trEnd-trStart)/fs;

        % find inclusions for trial start and end
        trStartIndx = find(spikeFrames >trStart,1, 'first');
        trEndIndx = find(spikeFrames <trEnd,1, 'last');

        % get the spikes
        trSpikes = spikeFrames(trStartIndx:  trEndIndx);
        trSpikes = trSpikes / fs; % convert to sec
        trSpikes = trSpikes- (stimOFF_Events{stimblk}(tr)/ fs); % rezero to alignment event

        % put into cell array for psth
        trSpikesStructOFF{tr} = trSpikes;
        trSpikesAll = [trSpikesAll; length(trSpikes)];
    end

    trSpikesMeanOff(stimblk) = sum(trSpikesAll)/length(stimOFF_Events{stimblk});

    %% PSTH

    trialLen            = mean(trialLenSec) * 1000;                % trial length ms
    edges               = 0:(binsize/1000): trialLenSec;
    binActual           = edges;
    nbins               = trialLen/binsize;                        % Bin duration in [ms]
    nobins              = 1000/binsize;                            % sec / bin
    all = [];
    ISI_off{stimblk} = [];
    trSpikesVarOffAll = [];
    for iTrial = 1:length(trSpikesStructOFF)
        all             = [all; trSpikesStructOFF{iTrial}];               % Concatenate spikes of all trials
        ISI_off{stimblk} = [ISI_off{stimblk}; diff(trSpikesStructOFF{iTrial})];
        trSpikesVarOffAll = [trSpikesVarOffAll; (trSpikesAll(iTrial) - trSpikesMeanOff(stimblk))^2];
        [trialHistTemp, ~] = histcounts(trSpikesStructOFF{iTrial},binActual);
        trialHistsoff{stimblk}(iTrial,:) = trialHistTemp;
    end

   
   ISIoff_mean(stimblk) = mean(ISI_off{stimblk});
   meanPSTH = mean(trialHistsoff{stimblk});
   OFFcountAverageSec     = (meanPSTH) * nobins;
   [A3(stimblk), peak_bin(stimblk)] = max(OFFcountAverageSec);
    A4(stimblk) = max(OFFcountAverageSec)*(1/exp(1));
    % Find the index after the peak where the PSTH drops below A4
    post_peak = OFFcountAverageSec(peak_bin(stimblk):end);
    below_A4 = find(post_peak < A4(stimblk), 1, 'first');
    if ~isempty(below_A4)
        tau_bin = peak_bin(stimblk) + below_A4 - 1; % Index in original array
        Ti_off(stimblk) = (tau_bin - peak_bin(stimblk)) * (binsize/1000); % Time difference (ms if binsize in ms)
    else
       below_A4 = OFFcountAverageSec(end);
       tau_bin = peak_bin(stimblk) + below_A4 - 1; % Index in original array
       Ti_off(stimblk) = (tau_bin - peak_bin(stimblk)) * (binsize/1000); % Time difference (ms if binsize in ms)
    end

    BI(stimblk) = (A1(stimblk)-A3(stimblk))/(A1(stimblk)+A3(stimblk));
    ratio(stimblk) = mean(ONcountAverageSec)/mean(OFFcountAverageSec);
    [~, pvalues(stimblk)] = ttest2(ONcountAverageSec,OFFcountAverageSec);
end

end