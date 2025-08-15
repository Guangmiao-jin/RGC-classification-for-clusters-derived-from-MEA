function data = processRetinaFlashStimDataFN1_V2(clusterFilepath)

%% defaults

%zScoreLower = 3; % minimum z score limit for responsivity

prestimTimeZScore = 1; % sec
stimTimeZScore = 1.5; % sec

prestimTimePSTH = 0.5; % sec
postStimTimePSTH = 1; % sec

responseQualityThreshold = 0.1; 

%% get the appropriate paths for stim on and off triggers
filepathPrefix = extractBefore(clusterFilepath, '_cluster');
stimOnFile = [filepathPrefix{:} '_triggerON.npy'];
stimOffFile = [filepathPrefix{:} '_triggerOFF.npy'];


%% load data
matSavePath = extractBefore(clusterFilepath, '.');
matSavePath = [matSavePath{:} '.mat'];

if ~exist(matSavePath)
    data = readHS2_FLAME(clusterFilepath);
    save(matSavePath, "data", "-v7.3");
else
    load(matSavePath);
end

%stimOnFrames = double(readNPY(stimOnFile));
%stimOffFrames = double(readNPY(stimOffFile));

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

%{
% % plotting for testing
blockOnStarts = [stimOnFrames(stimOnBreaks)];
blockOnEnds = [stimOnFrames(stimOnStopBreaks)];
scatter(stimOnFrames, ones(length(stimOnFrames))*100);
hold on
scatter(blockOnStarts, repmat(100.5, 1, length(blockOnStarts)));
scatter(blockOnEnds, repmat(100.5, 1, length(blockOnEnds)));
%}

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

%% make all the metrics we use to seperate out the cells
responseMetrics = createFlashMetrics(data, stimOnPerBlock, stimOffPerBlock, prestimTimeZScore, stimTimeZScore, prestimTimePSTH, postStimTimePSTH);

% find neuron ids contain having valid ISI
% and waveform
validNeurons = [];
for i = 1:length(responseMetrics.ISIon_mean)
 if  data.channelNames{6,i} >= 300
     validNeurons = [validNeurons; i];
 end
end

numValidClusters = length(validNeurons);
numNeurons = [];
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
        ISI_vio = responseMetrics.ISI_vio(i);
    if peaktimediff > 0 && ISI_vio < 0.05
        numNeurons = [numNeurons; i];
    end
    end
end



%respondingClustersIndex = max(responseMetrics.responseQuality,[],2)>responseQualityThreshold;
%respondingClusterIndxNum = find(respondingClustersIndex);

%% start the plotting

cellRasterFolder = extractBefore(clusterFilepath, '.');

binsize = 50; % ms

% set(0,'DefaultFigureVisible','off');
count = 0;

totalneurons.num_OnNeurons.trans = [];
totalneurons.num_OnNeurons.sus = [];

totalneurons.num_OffNeurons.trans = [];
totalneurons.num_OffNeurons.sus = [];

totalneurons.num_OnOffNeurons = [];
totalneurons.num_unconventional = [];

totalneurons.num_notRespond = [];

for i = 1:length(numNeurons)

    curCl = numNeurons(i);
    clusterID = curCl-1;

    count = count +1;

    disp(['On ' num2str(count) ' of ' num2str(length(numNeurons))]);

    clusterResponses.ZValueON = responseMetrics.ZValueON{curCl};
    clusterResponses.ZValueOFF = responseMetrics.ZValueOFF{curCl};
    clusterResponses.responseQuality = responseMetrics.responseQuality(curCl,:);
    clusterResponses.trialSpikes = responseMetrics.trialSpikes{curCl};
    clusterResponses.trialPSTHs = responseMetrics.trialPSTHs{curCl};
    clusterResponses.binEdges = responseMetrics.PSTH_binEdges;
    clusterResponses.BI = responseMetrics.BI{curCl};
    clusterResponses.TI_on =responseMetrics.TI_on{curCl};
    clusterResponses.TI_off = responseMetrics.TI_off{curCl};
    clusterResponses.ratio =responseMetrics.ratio{curCl};


   
    
        if max(clusterResponses.responseQuality)>responseQualityThreshold
            if mean((clusterResponses.BI)) >=0.33
                if mean(clusterResponses.TI_on,'omitnan') < 0.1
                    totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/trans';
                elseif mean(clusterResponses.TI_on,'omitnan') >= 0.1
                    totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/sus';
                end
            elseif mean((clusterResponses.BI)) <= - 0.33
                if mean(clusterResponses.TI_off,'omitnan') <= 0.1
                    totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/trans';
                elseif mean(clusterResponses.TI_off,'omitnan') >= 0.1
                    totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/sus';
                end
            elseif mean(clusterResponses.BI) < 0.33 && mean(clusterResponses.BI) > -0.33
                totalneurons.num_OnOffNeurons(end+1) = data.channelNames{4,curCl};
                cls = 'OnOffNeurons';
            else
                if mean(clusterResponses.ratio) >= 2
                    if mean(clusterResponses.TI_on,'omitnan') <= 0.1
                        totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                        cls = 'OnNeurons/trans';
                     elseif mean(clusterResponses.TI_on,'omitnan') >= 0.1
                        totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                        cls = 'OnNeurons/sus';
                    end
                elseif mean(clusterResponses.ratio) <= 0.5
                    if mean(clusterResponses.TI_off,'omitnan') <= 0.1
                        totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                        cls = 'OffNeurons/trans';
                    elseif mean(clusterResponses.TI_off,'omitnan') >= 0.1
                        totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                        cls = 'OffNeurons/sus';
                    end
                else
                    totalneurons.num_unconventional(end+1) = data.channelNames{4,curCl};
                    cls = 'unconventional';
                end
            end
        elseif max(clusterResponses.responseQuality)>0.05 && max(clusterResponses.responseQuality) <= responseQualityThreshold
            totalneurons.num_unconventional(end+1) = data.channelNames{4,curCl};
            cls = 'unconventional';
        else
            totalneurons.num_notRespond(end+1) = data.channelNames{4,curCl};
            cls = 'notRespond';
        end
   
        totalneurons.allneurons = [totalneurons.num_OnNeurons.trans, totalneurons.num_OnNeurons.sus, ...
            totalneurons.num_OffNeurons.trans, totalneurons.num_OffNeurons.sus,...
            totalneurons.num_OnOffNeurons, totalneurons.num_unconventional, totalneurons.num_notRespond];

     outputFolder = fullfile([cellRasterFolder{:} '_PSTHPlotsV2'], cls);
     if ~exist(outputFolder, 'dir')
         mkdir(outputFolder);
     end
    

    ax = plotAllRasterPSTHs_On_Off_ResponseV2(clusterResponses);

    sgtitle(['Cluster ID: ' num2str(clusterID) ' Spks: ' num2str(data.channelNames{6, curCl})]);
    tightfig;

    % saveName = sprintf('%s/cluster%04d.png',[cellRasterFolder{:} '_cellPlots'],i-1);
    saveName = sprintf('%s\\cluster%04d.png',outputFolder,clusterID);

    winHandle = gethwnd(gcf);
    cmndstr = sprintf('%s','MiniCap.exe -save ','"',saveName,'"',...
        ' -compress 9', ' -capturehwnd ', num2str(winHandle),' -exit');
    system(cmndstr);
    close
 %}
end

save(strcat(cellRasterFolder,"_responseMetrics.mat"),'responseMetrics');
save(strcat(cellRasterFolder,"_totalneuronsV2.mat"), 'totalneurons');  % 

% end
% set(0,'DefaultFigureVisible','on');
end
