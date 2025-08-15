function classifyNeuron(i, cellRasterFolder)

    figure('Units','normalized','Position',[0.1, 0.1, 0.8, 0.6]);
    subplot(1,2,1);
    % Load waveform image
    waveformPath = sprintf('%s/cluster%04d.png',[cellRasterFolder '_neurons_waveforms'],i-1);
    if exist(waveformPath, 'file')
        imshow(waveformPath);
    end

    % Load PTSH image
    subplot(1,2,2);
    ptshPath = sprintf('%s/cluster%04d.png',[cellRasterFolder '_PSTHPlots'],i-1);
    if exist(ptshPath, 'file')
        imshow(ptshPath);
    end

    drawnow;
    % User selection
end
