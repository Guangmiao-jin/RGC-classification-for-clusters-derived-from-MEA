function [ZvaluePerClusterBlkON, ZvaluePerClusterBlkOFF] = stimAlignedrank(spikeFrames, fs, stimON_Events, stimOFF_Events , postAlignTime)
binsize  =  50; %ms
preAlignTime = 1;
for stimblk = 1:length(stimON_Events)
    if stimblk < 4
    % alignment times in sec
    preAlignInFrames = preAlignTime * fs;
    postAlignInFrames = postAlignTime * fs;

    %% ON stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_on{stimblk} = [];
    for tr = 1:length(stimON_Events{stimblk})
        % get prestim and poststim spikes
        preStimStart = stimON_Events{stimblk}(tr)- preAlignInFrames;
        stimEnd = stimON_Events{stimblk}(tr)+ postAlignInFrames;

        trStartIndx = find(spikeFrames >preStimStart,1, 'first');
        trEndIndx = find(spikeFrames <stimEnd,1, 'last');

        stimSpikesON{stimblk}{tr} =spikeFrames(trStartIndx:  trEndIndx)-stimON_Events{stimblk}(tr);
        all_on{stimblk} = [all_on{stimblk} ; stimSpikesON{stimblk}{tr}];
    end
    trialLen            = (preAlignTime + postAlignTime) * 1000;   % trial length ms
    nbins               = trialLen/binsize;                        % Number of bins per trial 
    nobins              = 1000/binsize;                            % Number of bins per sec

   [counts, ~]     = histcounts(all_on{stimblk}/fs,nbins);
   ONcountAverageSec{stimblk}     = (counts/length(stimON_Events{stimblk})) * nobins;

    %% OFF stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_off{stimblk} = [];
    for tr = 1:length(stimOFF_Events{stimblk})
        % get prestim and poststim spikes
        preStimStart = stimOFF_Events{stimblk}(tr)- preAlignInFrames;
        stimEnd = stimOFF_Events{stimblk}(tr)+ postAlignInFrames;

        trStartIndx = find(spikeFrames >preStimStart,1, 'first');
        trEndIndx = find(spikeFrames <stimEnd,1, 'last');

        stimSpikesOFF{stimblk}{tr} =spikeFrames(trStartIndx:  trEndIndx)-stimOFF_Events{stimblk}(tr);
        all_off{stimblk} = [all_off{stimblk} ; stimSpikesOFF{stimblk}{tr}];
    end
    trialLen            = (preAlignTime + postAlignTime) * 1000;   % trial length ms
    nbins               = trialLen/binsize;                        % Number of bins per trial 
    nobins              = 1000/binsize;                            % Number of bins per sec

    [counts, ~]     = histcounts(all_off{stimblk}/fs,nbins);
    OFFcountAverageSec{stimblk}     = (counts/length(stimOFF_Events{stimblk})) * nobins;    
 
    %% no stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if stimblk == 1
            nostimSpikes{stimblk} = spikeFrames(1 < spikeFrames & spikeFrames < 2*60*fs)/fs;
            trialLen            = 2 * 1000;            % trial length ms
            nbins               = trialLen/binsize;                        % Number of bins per trial 
            nobins              = 1000/binsize;                            % Number of bins per sec
             [counts, ~]     = histcounts(nostimSpikes{stimblk}/fs,nbins);
             NoStimcountAverageSec{stimblk}     = (counts/120) * nobins;
        elseif stimblk < 4
            nostimSpikes{stimblk} = (spikeFrames(stimOFF_Events{stimblk-1}(end) < spikeFrames & spikeFrames < stimOFF_Events{stimblk-1}(end) + 2*60*fs )-stimOFF_Events{stimblk-1}(end))/fs;
            trialLen            = 2 * 1000;            % trial length ms
            nbins               = trialLen/binsize;                        % Number of bins per trial 
            nobins              = 1000/binsize;                            % Number of bins per sec
            [counts, ~]     = histcounts(nostimSpikes{stimblk}/fs,nbins);
            NoStimcountAverageSec{stimblk}     = (counts/120) * nobins;
        end
     ZvaluePerClusterBlkON(stimblk)= (mean(ONcountAverageSec{stimblk}(21:40))-mean(ONcountAverageSec{stimblk}(1:20)))/ std(ONcountAverageSec{stimblk}(1:20));
     ZvaluePerClusterBlkOFF(stimblk) = (mean(OFFcountAverageSec{stimblk}(21:40))-mean(OFFcountAverageSec{stimblk}(1:20)))/ std(OFFcountAverageSec{stimblk}(1:20));
     %[offpval(stimblk),~,~] = kruskalwallis(NoStimcountAverageSec{stimblk},OFFcountAverageSec{stimblk});
     %pval(stimblk) = ;


    end
end

end