function [A1, A2, BI, ISI_on, ISI_off, Ti_on, Ti_off] = find_On_Off_Response(spikeFrames, binsize, fs, stimON_Events, stimOFF_Events, preAlignTime, postAlignTime)
%% ON stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stimblk = 1:length(stimON_Events)
    % alignment times in sec
    preAlignInFrames = preAlignTime * fs;
    postAlignInFrames = postAlignTime * fs;

    % for each trial
    for tr = 1:length(stimON_Events{stimblk})

        % trial times in frames
        trStart = stimON_Events{stimblk}(tr)- preAlignInFrames;
        trEnd = stimON_Events{stimblk}(tr)+ postAlignInFrames;

        % find inclusions for trial start and end
        trStartIndx = find(spikeFrames >trStart,1, 'first');
        trEndIndx = find(spikeFrames <trEnd,1, 'last');

        % get the spikes
        trSpikes = spikeFrames(trStartIndx:  trEndIndx);
        trSpikes = trSpikes / fs; % convert to sec
        trSpikes = trSpikes- (stimON_Events{stimblk}(tr)/ fs); % rezero to alignment event

        % put into cell array for psth
        trSpikesStructON{tr} = trSpikes;
    end

    %% PSTH

    all = [];
    ISI_on{stimblk} = [];
    for iTrial = 1:length(trSpikesStructON)
        all             = [all; trSpikesStructON{iTrial}];               % Concatenate spikes of all trials
        ISI_on{stimblk} = [ISI_on{stimblk}; diff(trSpikesStructON{iTrial})];
    end

    trialLen            = (preAlignTime + postAlignTime) * 1000;    % trial length ms
    nbins               = trialLen/binsize;                        % Bin duration in [ms]
    nobins              = 1000/binsize;                            % No of bins/sec

    [counts, ~]     = histcounts(all,nbins);
    ONcountAverageSec     = (counts/length(stimON_Events{stimblk})) * nobins;
    A1(stimblk) = max(ONcountAverageSec);
    steady_rate(stimblk) = mean(ONcountAverageSec(end-25:end));
    Ti_on(stimblk) = (A1(stimblk)-steady_rate(stimblk))/A1(stimblk);
end

%% OFF stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stimblk = 1:length(stimOFF_Events)
    % alignment times in sec
    preAlignInFrames = preAlignTime * fs;
    postAlignInFrames = postAlignTime * fs;

    % for each trial
    for tr = 1:length(stimOFF_Events{stimblk})

        % trial times in frames
        trStart = stimOFF_Events{stimblk}(tr)- preAlignInFrames;
        trEnd = stimOFF_Events{stimblk}(tr)+ postAlignInFrames;

        % find inclusions for trial start and end
        trStartIndx = find(spikeFrames >trStart,1, 'first');
        trEndIndx = find(spikeFrames <trEnd,1, 'last');

        % get the spikes
        trSpikes = spikeFrames(trStartIndx:  trEndIndx);
        trSpikes = trSpikes / fs; % convert to sec
        trSpikes = trSpikes- (stimOFF_Events{stimblk}(tr)/ fs); % rezero to alignment event

        % put into cell array for psth
        trSpikesStructOFF{tr} = trSpikes;
    end
    %% PSTH

    all = [];
    ISI_off{stimblk} = [];
    for iTrial = 1:length(trSpikesStructOFF)
        all             = [all; trSpikesStructOFF{iTrial}];               % Concatenate spikes of all trials
        ISI_off{stimblk}= [ISI_off{stimblk}; diff(trSpikesStructOFF{iTrial})];
    end

    trialLen            = (preAlignTime + postAlignTime) * 1000;    % trial length ms
    nbins               = trialLen/binsize;                        % Bin duration in [ms]
    nobins              = 1000/binsize;                            % No of bins/sec

    [counts, ~]     = histcounts(all,nbins);
    OFFcountAverageSec     = (counts/length(stimOFF_Events{stimblk})) * nobins;
    A2(stimblk) = max(OFFcountAverageSec);
    steady_rate(stimblk) = mean(OFFcountAverageSec(end-25:end));
    Ti_off(stimblk) = (A2(stimblk)-steady_rate(stimblk))/A2(stimblk);

    BI(stimblk) = (A1(stimblk)-A2(stimblk))/(A1(stimblk)+A2(stimblk));
end
end