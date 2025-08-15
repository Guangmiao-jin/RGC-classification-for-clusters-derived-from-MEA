function processRetinaFlashStimDataFN1conditions(clusterFilepath)
%% -------- load data --------
responseQualityThreshold = 0.1;

preSavePath = extractBefore(clusterFilepath, '.');
matSavePath = [preSavePath{:} '.mat'];
load(matSavePath,'data');

responsePath = [preSavePath{:} '_responseMetrics.mat'];
load(responsePath,'responseMetrics');

%% -------- select valid neurons --------
validNeurons = [];
for i = 1:length(responseMetrics.ISIon_mean)
    if  data.channelNames{6,i} >= 300
        validNeurons = [validNeurons; i]; %#ok<AGROW>
    end
end

numValidClusters = length(validNeurons);
numNeurons = [];
for idx = 1:numValidClusters
    i = validNeurons(idx);
    meanwaveforms = data.waveformClusterMeans(i,:);
    
    % positive and negative peaks
    [invertedPeaks, t0] = findpeaks(-meanwaveforms);
    [~,t1] = max(meanwaveforms);

    if ~isempty(invertedPeaks)
        numNegativeVally = sum(invertedPeaks>0);
        if numNegativeVally == 1 && abs(invertedPeaks(1)) > 1000
            peaktimediff = t1 - t0(1);
            ISI_vio = responseMetrics.ISI_vio(i);
            if peaktimediff > 0 && ISI_vio < 0.05
                numNeurons = [numNeurons; i];
            end
        end
    end
end


%% -------- classification --------
condition_names = {'scotopic','mesopic','photopic','all'};
results = struct([]);

count = 0;

for k = 1:numel(condition_names)
    condName = condition_names{k};

    % 全部初始化为数值矩阵（列向量）
    N = struct();
    N.name = condName;  % 只是个标签，不参与矩阵存储
    N.num_OnNeurons = struct('trans', [], 'sus', []);
    N.num_OffNeurons = struct('trans', [], 'sus', []);
    N.num_OnOffNeurons = [];
    N.num_unconventional = [];
    N.num_notRespond = [];

    for i = 1:numel(numNeurons)
        curCl = numNeurons(i);

        count = count + 1;
        disp(['On ' num2str(count) ' of ' num2str(4*numel(numNeurons))]);

        % 只取分类所需的度量
        clusterResponses.responseQuality  = responseMetrics.responseQuality(curCl,:);
        clusterResponses.BI               = responseMetrics.BI{curCl};
        clusterResponses.TI_on            = responseMetrics.TI_on{curCl};
        clusterResponses.TI_off           = responseMetrics.TI_off{curCl};
        clusterResponses.ratio            = responseMetrics.ratio{curCl};

        rq_max = max(clusterResponses.responseQuality);

        if k ~= 4  % 单条件
            BI_k    = clusterResponses.BI(k);
            ratio_k = clusterResponses.ratio(k);
        else       % "all"：跨条件均值
            BI_k    = mean(clusterResponses.BI,   'omitnan');
            ratio_k = mean(clusterResponses.ratio,'omitnan');
        end

        TIon  = mean(clusterResponses.TI_on,  'omitnan');
        TIoff = mean(clusterResponses.TI_off, 'omitnan');

        % ⚠️ 假设这里拿到的是**数值**ID（如通道号）；若不是数值，请改成数值来源
        ch = data.channelNames{4, curCl};

        % 分类逻辑（保持你的阈值）
        if rq_max > responseQualityThreshold

            if BI_k >= 0.33
                if TIon < 0.1
                    N.num_OnNeurons.trans(end+1,1) = ch;
                else
                    N.num_OnNeurons.sus(end+1,1)   = ch;
                end

            elseif BI_k <= -0.33
                if TIoff <= 0.1
                    N.num_OffNeurons.trans(end+1,1) = ch;
                else
                    N.num_OffNeurons.sus(end+1,1)   = ch;
                end

            elseif (BI_k > -0.33) && (BI_k < 0.33)
                N.num_OnOffNeurons(end+1,1) = ch;

            else
                % ratio 兜底
                if ratio_k >= 2
                    if TIon <= 0.1
                        N.num_OnNeurons.trans(end+1,1) = ch;
                    else
                        N.num_OnNeurons.sus(end+1,1)   = ch;
                    end
                elseif ratio_k <= 0.5
                    if TIoff <= 0.1
                        N.num_OffNeurons.trans(end+1,1) = ch;
                    else
                        N.num_OffNeurons.sus(end+1,1)   = ch;
                    end
                else
                    N.num_unconventional(end+1,1) = ch;
                end
            end

        elseif (rq_max > 0.05) && (rq_max <= responseQualityThreshold)
            N.num_unconventional(end+1,1) = ch;
        else
            N.num_notRespond(end+1,1) = ch;
        end
    end

    % 汇总到一个列向量（数值）
    N.S = [ N.num_OnNeurons.trans; N.num_OnNeurons.sus; ...
            N.num_OffNeurons.trans; N.num_OffNeurons.sus; ...
            N.num_OnOffNeurons; N.num_unconventional; N.num_notRespond ];

    % 计数
    N.counts = [ numel(N.num_OnNeurons.trans); ...
                 numel(N.num_OnNeurons.sus); ...
                 numel(N.num_OffNeurons.trans); ...
                 numel(N.num_OffNeurons.sus); ...
                 numel(N.num_OnOffNeurons); ...
                 numel(N.num_unconventional); ...
                 numel(N.num_notRespond) ];

    results{k} = N; 
end

% 保存结果（包含四个条件）
save(strcat(preSavePath, "_totalneuronsV3.mat"), 'results');

% ------------ Build sets (A,B,C,D) with OnOff included, robustly ------------
% ------------ Build sets (A,B,C,D) with OnOff included, robustly ------------
onSets  = cell(1,4);
offSets = cell(1,4);

% 工具：把任意数值向量 -> 行向量集合（去重、去 NaN）
mkset = @(v) unique(v(~isnan(v)));
row   = @(v) v(:).';

for k = 1:4
    % ON = On_trans ∪ On_sus   （不计 OnOff）
    onSets{k}  = row( mkset([results{k}.num_OnNeurons.trans; ...
                              results{k}.num_OnNeurons.sus]) );

    % OFF = Off_trans ∪ Off_sus （不计 OnOff）
    offSets{k} = row( mkset([results{k}.num_OffNeurons.trans; ...
                              results{k}.num_OffNeurons.sus]) );
end

condLabels = {'scotopic','mesopic','photopic','all'};

% ---- 自检：集合大小（来自集合本身）----
onTotals_sets  = cellfun(@numel, onSets);
offTotals_sets = cellfun(@numel, offSets);
fprintf('\n[SET SIZES from sets]\n');
for k = 1:4
    fprintf(' %s: ON=%d, OFF=%d\n', condLabels{k}, onTotals_sets(k), offTotals_sets(k));
end

% ---- 16 区域计数（稳健矢量化；不会出现全 0）----
onCounts16  = membership4count(onSets{:});
offCounts16 = membership4count(offSets{:});

% ---- 边缘总数（由 16 区域回推），应与上面一致 ----
onTotals_fromCounts  = totalsFromCounts(onCounts16);
offTotals_fromCounts = totalsFromCounts(offCounts16);

fprintf('\n[SET SIZES from 16-region counts]\n');
fprintf(' ON [A B C D] = [%d %d %d %d]\n',  onTotals_fromCounts);
fprintf('OFF [A B C D] = [%d %d %d %d]\n',  offTotals_fromCounts);

% 一致性检查（若不一致会警告）
checkTotals('ON',  condLabels, onTotals_sets,  onTotals_fromCounts);
checkTotals('OFF', condLabels, offTotals_sets, offTotals_fromCounts);

% ---- 可选：验证“所有区域之和 = 并集大小” ----
assert(sum(onCounts16)  == numel(unique([onSets{1}, onSets{2}, onSets{3}, onSets{4}])),  'ON 16-region sum != union size');
assert(sum(offCounts16) == numel(unique([offSets{1}, offSets{2}, offSets{3}, offSets{4}])), 'OFF 16-region sum != union size');

% ---- 画图（保留你的函数样式，圆旁标总数）----
figure('Name','ON 4-set Venn');
plotVenn4_withTotals(onCounts16, condLabels, 'ON overlap (no OnOff)');

figure('Name','OFF 4-set Venn');
plotVenn4_withTotals(offCounts16, condLabels, 'OFF overlap (no OnOff)');

% ===================== helpers =====================
function counts16 = membership4count(A,B,C,D)
% A,B,C,D: 行/列向量数值集合；按 ABCD 生成 0..15 掩码并累计到 16 区域。
    A = A(:); B = B(:); C = C(:); D = D(:);
    U = unique([A;B;C;D]);           % 宇集（至少包含所有进集合的元素）
    if isempty(U)
        counts16 = zeros(16,1);
        return;
    end
    inA = ismember(U, A);
    inB = ismember(U, B);
    inC = ismember(U, C);
    inD = ismember(U, D);
    % 掩码：bit1=A, bit2=B, bit3=C, bit4=D
    code = uint8(inA) + 2*uint8(inB) + 4*uint8(inC) + 8*uint8(inD); % 0..15
    counts16 = accumarray(double(code)+1, 1, [16 1]);               % 索引=code+1
end

function totals = totalsFromCounts(counts16)
% 由 16 区域计数求边缘总数 |A|,|B|,|C|,|D|；不做 +1。
    totals = zeros(1,4);
    for bit = 1:4
        idx = find(bitand(0:15, 2^(bit-1)) > 0); % 哪些区域包含该集合
        totals(bit) = sum(counts16(idx));
    end
end

function checkTotals(tag, labels, totals_sets, totals_counts)
    mismatch = totals_sets(:) ~= totals_counts(:);
    if any(mismatch)
        for j = find(mismatch).'
            warning('%s "%s": totals mismatch — sets=%d, fromCounts=%d', ...
                tag, labels{j}, totals_sets(j), totals_counts(j));
        end
    end
end

function plotVenn4_withTotals(counts16, labels, titleStr)
    if numel(labels)~=4, labels = {'A','B','C','D'}; end

    cla reset; ax = gca; hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

    % 四个圆的布局（A左、B右、C上、D下）
    r = 2.2;
    ctr = [ -1.8  0.0;   % A
             1.8  0.0;   % B
             0.0  1.8;   % C
             0.0 -1.8 ]; % D
    th = linspace(0,2*pi,400);
    for i = 1:4
        plot(ax, ctr(i,1)+r*cos(th), ctr(i,2)+r*sin(th), 'LineWidth',1.5);
    end

    % ---- 动态计算各区域位置 ----
    Uc   = mean(ctr,1);
    pos  = zeros(16,2);
    for code = 0:15
        S = find(bitget(code,1:4));
        if isempty(S)               % 0000 区域：我们不想显示
            pos(code+1,:) = [NaN NaN]; 
        else
            cS = mean(ctr(S,:),1);
            k  = numel(S);
            if k==1                 % 单集：往外推
                d = cS - Uc; if norm(d)==0, d=[1 0]; end
                pos(code+1,:) = cS + 0.65*r*(d/norm(d));
            elseif k==2             % 两两交
                pos(code+1,:) = mean(ctr(S,:),1);
            elseif k==3             % 三重交
                pos(code+1,:) = 0.85*mean(ctr(S,:),1) + 0.15*Uc;
            else                    % 四重交
                pos(code+1,:) = Uc + [0, 0.2];
            end
        end
    end

    % ---- 画区域数字（跳过 code==0）----
    for code = 1:15
        text(ax, pos(code+1,1), pos(code+1,2), num2str(counts16(code+1)), ...
            'HorizontalAlignment','center','FontSize',10);
    end

    % ---- 计算边缘总数并在圆外显示 ----
    totals = zeros(1,4);
    for bit = 1:4
        idx = find(bitand(0:15, 2^(bit-1))>0);
        totals(bit) = sum(counts16(idx));
    end
    for i = 1:4
        d = ctr(i,:) - Uc; if norm(d)==0, d=[1 0]; end
        p_label = ctr(i,:) + (r*1.35)*(d/norm(d));
        text(ax, p_label(1), p_label(2), ...
            sprintf('%s (n=%d)', labels{i}, totals(i)), ...
            'FontWeight','bold','HorizontalAlignment','center', ...
            'Margin',3,'Clipping','on');
    end

    % 固定视野
    xlim(ax, [min(ctr(:,1))-1.7*r, max(ctr(:,1))+1.7*r]);
    ylim(ax, [min(ctr(:,2))-1.7*r, max(ctr(:,2))+1.7*r]);

    title(ax, sprintf('%s — union=%d (sum of 15 regions)', titleStr, sum(counts16(2:end))));
    hold(ax,'off');
end


end