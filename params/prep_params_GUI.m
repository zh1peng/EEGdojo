function prep_params_GUI
    % 创建主 UI 窗口
    fig = uifigure('Name', 'EEGdojo🥋 Parameters Setup', ...
                  'Position', [200 200 1300 700], ...
                   'Color', [0.65 0.76 1]); % 浅蓝色背景
    
    uilabel(fig, 'Text', 'EEGDojo🥋', ...
        'FontWeight', 'bold', ...
        'FontSize', 24, ... % 根据需要调整字体大小
        'Position', [550 660 400 40]);
    % "Clickable Link" as a button
uihyperlink(fig, 'Text', 'Github: https://github.com/zh1peng/EEGdojo', ...
    'URL', 'https://github.com/zh1peng/EEGdojo', ...
    'Position', [700 655 200 30], ...
    'FontSize', 12);

    % 创建滚动容器
    scrollPanel = uipanel(fig, 'Position', [0 10 1280 650], ...
                          'Scrollable', 'on', ...
                          'BackgroundColor', [0.65 0.76 1]);

    % 依次创建各个模块
    createChannelInfoPanel(scrollPanel);
    createDownsamplingPanel(scrollPanel);
    createEpochingPanel(scrollPanel);
    createBadChannelPanel(scrollPanel);
    createReReferencingPanel(scrollPanel);
    createICAParametersPanel(scrollPanel);
    createBadPointDetectionPanel(scrollPanel);
    createEpochRestPanel(scrollPanel);

    % 创建 "Generate Params" 按钮
    generateBtn = uibutton(scrollPanel, 'push', ...
        'Text', 'Generate Params', ...
        'Position', [850 120 200 80], ...
        'BackgroundColor', [0.2 0.6 1], ... % 蓝色按钮背景
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(btn,event) generateParamsCallback(fig));

    %% 创建各个面板的函数

    function createChannelInfoPanel(parent)
        % 创建 "Channel Information" 面板
        chanPanel = uipanel(parent, 'Title', 'Channel Information', ...
            'Position', [10 800 1250 180], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % "All Channels" 输入框
        uilabel(chanPanel, 'Text', 'All Channels (e.g., 1:129):', ...
            'Position', [10 130 180 22], 'FontSize', 11);
        chanAllEdit = uieditfield(chanPanel, 'text', ...
            'Position', [200 130 200 22], ...
            'Value', '1:129', ...
            'Tag', 'chanAllEdit');

        % "Reference Channel" 输入框
        uilabel(chanPanel, 'Text', 'Reference Channel:', ...
            'Position', [10 90 180 22], 'FontSize', 11);
        refChanEdit = uieditfield(chanPanel, 'numeric', ...
            'Position', [200 90 200 22], ...
            'Value', 129, ...
            'Tag', 'refChanEdit');

        % "EOG Channel Index" 输入框
        uilabel(chanPanel, 'Text', 'EOG Channel Index:', ...
            'Position', [10 50 180 22], 'FontSize', 11);
        eogIdxEdit = uieditfield(chanPanel, 'text', ...
            'Position', [200 50 200 22], ...
            'Value', '[]', ...
            'Tag', 'eogIdxEdit');

        % "EOG Channel Labels" 输入框
        uilabel(chanPanel, 'Text', 'EOG Channel Labels:', ...
            'Position', [10 10 180 22], 'FontSize', 11);
        eogLabelsEdit = uieditfield(chanPanel, 'text', ...
            'Position', [200 10 200 22], ...
            'Value', '[]', ...
            'Tag', 'eogLabelsEdit');

        % "Known Bad Channel Index" 输入框
        uilabel(chanPanel, 'Text', 'Known Bad Channel Index:', ...
            'Position', [420 130 200 22], 'FontSize', 11);
        knownBadIdxEdit = uieditfield(chanPanel, 'text', ...
            'Position', [630 130 200 22], ...
            'Value', '[]', ...
            'Tag', 'knownBadIdxEdit');

        % "Known Bad Channel Labels" 输入框
        uilabel(chanPanel, 'Text', 'Known Bad Channel Labels:', ...
            'Position', [420 90 200 22], 'FontSize', 11);
        knownBadLabelsEdit = uieditfield(chanPanel, 'text', ...
            'Position', [630 90 200 22], ...
            'Value', '[]', ...
            'Tag', 'knownBadLabelsEdit');

        % Additional instructions
        % uilabel(chanPanel, 'Text', 'Specify channel indices and labels as vectors (e.g., [1 2 3]).', ...
        %     'Position', [10 160 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createDownsamplingPanel(parent)
        % 创建 "Downsampling and Filtering" 面板
        downPanel = uipanel(parent, 'Title', 'Downsampling and Filtering Parameters', ...
            'Position', [10 610 620 180], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Downsampling Rate
        uilabel(downPanel, 'Text', 'Downsampling Rate (Hz):', ...
            'Position', [10 130 180 22], 'FontSize', 11);
        downRateEdit = uieditfield(downPanel, 'numeric', ...
            'Position', [200 130 100 22], ...
            'Value', 250, ...
            'Tag', 'downRateEdit');

        % Low Cutoff
        uilabel(downPanel, 'Text', 'Low Cutoff (Hz):', ...
            'Position', [10 90 150 22], 'FontSize', 11);
        lowCutEdit = uieditfield(downPanel, 'numeric', ...
            'Position', [200 90 100 22], ...
            'Value', 0.1, ...
            'Tag', 'lowCutEdit');

        % High Cutoff
        uilabel(downPanel, 'Text', 'High Cutoff (Hz):', ...
            'Position', [10 50 150 22], 'FontSize', 11);
        highCutEdit = uieditfield(downPanel, 'numeric', ...
            'Position', [200 50 100 22], ...
            'Value', 40, ...
            'Tag', 'highCutEdit');

        % Power Line Frequency
        uilabel(downPanel, 'Text', 'Power Line Frequency (Hz):', ...
            'Position', [10 10 180 22], 'FontSize', 11);
        plfEdit = uieditfield(downPanel, 'numeric', ...
            'Position', [200 10 100 22], ...
            'Value', 60, ...
            'Tag', 'plfEdit');

        % Additional instructions
        % uilabel(downPanel, 'Text', 'Enter numerical values for filtering parameters.', ...
        %     'Position', [320 130 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createEpochingPanel(parent)
        % 创建 "Epoching Parameters" 面板
        epochPanel = uipanel(parent, 'Title', 'Epoching Parameters', ...
           'Position', [640 610 620 180], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Epoch Task Marker
        uilabel(epochPanel, 'Text', 'Epoch Task Marker:', ...
            'Position', [10 90 180 22], 'FontSize', 11);
        markerEdit = uieditfield(epochPanel, 'text', ...
            'Position', [200 90 200 22], ...
            'Value', 'instructed_toCloseEyes', ...
            'Tag', 'markerEdit');

        % Epoch Window
        uilabel(epochPanel, 'Text', 'Epoch Window [start, end] (s):', ...
            'Position', [10 50 200 22], 'FontSize', 11);
        epochWindowEdit = uieditfield(epochPanel, 'text', ...
            'Position', [220 50 150 22], ...
            'Value', '[1, 39]', ...
            'Tag', 'epochWindowEdit');

        % Additional instructions
        % uilabel(epochPanel, 'Text', 'Specify window in seconds as a vector (e.g., [1 39]).', ...
        %     'Position', [10 120 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createBadChannelPanel(parent)
        % 创建 "Bad Channel Detection" 面板
        badChanPanel = uipanel(parent, 'Title', 'Bad Channel Detection Parameters', ...
            'Position', [10 220 620 380], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Reject Channel Normalization
        rejectNormChk = uicheckbox(badChanPanel, 'Text', 'Reject Channel Normalization', ...
            'Position', [10 330 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'rejectNormChk');

        % Specific Frequency Range
        uilabel(badChanPanel, 'Text', 'Specific Frequency Range (Hz):', ...
            'Position', [270 330 200 22], 'FontSize', 11);
        specFreqEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [530 330 80 22], ...
            'Value', '[0.1 40]', ...
            'Tag', 'specFreqEdit');

        % Specific Threshold
        uilabel(badChanPanel, 'Text', 'Specific Threshold (z-score):', ...
            'Position', [270 300 200 22], 'FontSize', 11);
        specThreshEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [530 300 80 22], ...
            'Value', '[-5 5]', ...
            'Tag', 'specThreshEdit');

        % Enable Specific Frequency Range Checkbox
        specFreqChk = uicheckbox(badChanPanel, 'Text', 'Enable Specific Frequency Range', ...
            'Position', [10 300 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'specFreqChk');

        % Enable Specific Threshold Checkbox
        specThreshChk = uicheckbox(badChanPanel, 'Text', 'Enable Specific Threshold', ...
            'Position', [10 270 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'specThreshChk');

        % Reject Channel Threshold
        uilabel(badChanPanel, 'Text', 'Reject Channel Threshold (z-score):', ...
            'Position', [270 270 250 22], 'FontSize', 11);
        rejectThreshEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [530 270 80 22], ...
            'Value', '[-5 5]', ...
            'Tag', 'rejectThreshEdit');

        % Kurtosis Based Rejection
        kurtOnChk = uicheckbox(badChanPanel, 'Text', 'Enable Kurtosis Based Rejection', ...
            'Position', [10 240 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'kurtOnChk');
        uilabel(badChanPanel, 'Text', 'Kurtosis Threshold (z-score):', ...
            'Position', [270 240 200 22], 'FontSize', 11);
        kurtThreshEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [530 240 80 22], ...
            'Value', '[-5 5]', ...
            'Tag', 'kurtThreshEdit');

        % Probability Based Rejection
        probOnChk = uicheckbox(badChanPanel, 'Text', 'Enable Probability Based Rejection', ...
            'Position', [10 210 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'probOnChk');
        uilabel(badChanPanel, 'Text', 'Probability Threshold (z-score):', ...
            'Position', [270 210 200 22], 'FontSize', 11);
        probThreshEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [530 210 80 22], ...
            'Value', '[-5 5]', ...
            'Tag', 'probThreshEdit');

        % FASTER Channel Variance
        fasterVarChk = uicheckbox(badChanPanel, 'Text', 'Enable FASTER Channel Variance', ...
            'Position', [10 180 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'fasterVarChk');
        uilabel(badChanPanel, 'Text', 'FASTER Channel Variance Threshold:', ...
            'Position', [270 180 250 22], 'FontSize', 11);
        fasterVarThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 180 80 22], ...
            'Value', 3, ...
            'Tag', 'fasterVarThreshEdit');

        % FASTER Channel Mean Correlation
        fasterMeanCorrChk = uicheckbox(badChanPanel, 'Text', 'Enable FASTER Channel Mean Correlation', ...
            'Position', [10 150 300 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'fasterMeanCorrChk');
        uilabel(badChanPanel, 'Text', 'FASTER Channel Mean Correlation Threshold:', ...
            'Position', [270 150 300 22], 'FontSize', 11);
        fasterMeanCorrThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 150 80 22], ...
            'Value', 3, ...
            'Tag', 'fasterMeanCorrThreshEdit');

        % FASTER Channel Hurst
        fasterHurstChk = uicheckbox(badChanPanel, 'Text', 'Enable FASTER Channel Hurst', ...
            'Position', [10 120 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'fasterHurstChk');
        uilabel(badChanPanel, 'Text', 'FASTER Channel Hurst Threshold:', ...
            'Position', [270 120 250 22], 'FontSize', 11);
        fasterHurstThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 120 80 22], ...
            'Value', 3, ...
            'Tag', 'fasterHurstThreshEdit');

        % % Clean Raw Data Parameters
        % uilabel(badChanPanel, 'Text', 'Clean Raw Data Parameters:', ...
        %     'Position', [650 220 200 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);


         % Clean Raw Data Clean Drifts High Pass
        uilabel(badChanPanel, 'Text', 'Clean Drifts High Pass [low, high] (Hz):', ...
            'Position', [10 90 300 22], 'FontSize', 11);
        cleanDriftsHighPassEdit = uieditfield(badChanPanel, 'text', ...
            'Position', [270 90 160 22], ...
            'Value', '[0.25, 0.75]', ...
            'Tag', 'cleanDriftsHighPassEdit');

        % Clean Raw Data Flat Line On
        cleanFlatLineChk = uicheckbox(badChanPanel, 'Text', 'Enable Flat Line Detection', ...
            'Position', [10 60 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'cleanFlatLineChk');
        uilabel(badChanPanel, 'Text', 'Flat Line Threshold (seconds):', ...
            'Position', [270 60 200 22], 'FontSize', 11);
        cleanFlatLineThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 60 80 22], ...
            'Value', 5, ...
            'Tag', 'cleanFlatLineThreshEdit');

       

        % Clean Raw Data Clean Channel
        cleanChanChk = uicheckbox(badChanPanel, 'Text', 'Enable Channel Cleaning', ...
            'Position', [10 30 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'cleanChanChk');

        % Clean Raw Data Line Noise Threshold
        uilabel(badChanPanel, 'Text', 'Line Noise Threshold:', ...
            'Position', [270 30 150 22], 'FontSize', 11);
        cleanLineNoiseThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 30 80 22], ...
            'Value', 4, ...
            'Tag', 'cleanLineNoiseThreshEdit');

        % Clean Raw Data Correlation Threshold
        uilabel(badChanPanel, 'Text', 'Correlation Threshold:', ...
            'Position', [270 10 150 22], 'FontSize', 11);
        cleanCorrThreshEdit = uieditfield(badChanPanel, 'numeric', ...
            'Position', [530 10 80 22], ...
            'Value', 0.6, ...
            'Tag', 'cleanCorrThreshEdit');

        % % Additional instructions
        % uilabel(badChanPanel, 'Text', 'Configure bad channel detection criteria.', ...
        %     'Position', [10 250 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createReReferencingPanel(parent)
        % 创建 "Re-referencing Parameters" 面板
        reRefPanel = uipanel(parent, 'Title', 'Re-referencing Parameters', ...
            'Position', [640 520 620 80], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Channel to Average Checkbox
        chanAvgChk = uicheckbox(reRefPanel, 'Text', 'Enable Channel to Average Re-referencing', ...
            'Position', [10 30 300 22], 'FontSize', 11, 'Value', false, ...
            'Tag', 'chanAvgChk');

        % Channels to Average
        uilabel(reRefPanel, 'Text', 'Channels to Average:', ...
            'Position', [10 10 150 22], 'FontSize', 11);
        chanAvgEdit = uieditfield(reRefPanel, 'text', ...
            'Position', [270 10 200 22], ...
            'Value', '[]', ...
            'Tag', 'chanAvgEdit');

        % % Additional instructions
        % uilabel(reRefPanel, 'Text', 'Specify channels as vectors (e.g., [1 2 3]).', ...
        %     'Position', [700 100 300 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createICAParametersPanel(parent)
        % 创建 "ICA Parameters" 面板
        icaPanel = uipanel(parent, 'Title', 'ICA Parameters', ...
            'Position', [640 330 620 180], ...
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Filter ICA On
        filterICAonChk = uicheckbox(icaPanel, 'Text', 'Filter ICA On', ...
            'Position', [10 130 150 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'filterICAonChk');
           % ICA Low Cutoff
        uilabel(icaPanel, 'Text', 'ICA Low Cutoff (Hz):', ...
            'Position',  [130 130 120 22], 'FontSize', 11);
        icaLowCutEdit = uieditfield(icaPanel, 'numeric', ...
            'Position',[270 130 100 22], ...
            'Value', 1, ...
            'Tag', 'icaLowCutEdit');

         
        % IC Label On
        icLabelChk = uicheckbox(icaPanel, 'Text', 'IC Label On', ...
            'Position', [10 100 150 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'icLabelChk');
      
         % IC Label Threshold
        uilabel(icaPanel, 'Text', 'IC Label Threshold:', ...
            'Position', [130 100 150 22], 'FontSize', 11);
        icLabelThreshEdit = uieditfield(icaPanel, 'text', ...
            'Position', [270 100 300 22], ...
            'Value', '[NaN NaN; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; NaN NaN]', ...
            'Tag', 'icLabelThreshEdit');

        % FASTER On
        fasterChk = uicheckbox(icaPanel, 'Text', 'FASTER On', ...
            'Position', [10 70 120 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'fasterChk');

       
% Number of Runs
        uilabel(icaPanel, 'Text', 'Number of Runs:', ...
            'Position', [10 40 120 22], 'FontSize', 11);
        nRunEdit = uieditfield(icaPanel, 'numeric', ...
            'Position', [130 40 100 22], ...
            'Value', 1, ...
            'Tag', 'nRunEdit');
       

        % % Additional instructions
        % uilabel(icaPanel, 'Text', 'Configure ICA processing parameters.', ...
        %     'Position', [10 250 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    function createBadPointDetectionPanel(parent)
        % 创建 "Bad Point Detection Parameters" 面板
        badPointPanel = uipanel(parent, 'Title', 'Bad Point Detection Parameters', ...
            'Position', [640 220 620 100], ... % 调整位置以适应滚动
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Bad Point Detection On
        badPointOnChk = uicheckbox(badPointPanel, 'Text', 'Enable Bad Point Detection', ...
            'Position', [10 60 250 22], 'FontSize', 11, 'Value', true, ...
            'Tag', 'badPointOnChk');

        % Burst Criterion
        uilabel(badPointPanel, 'Text', 'Burst Criterion:', ...
            'Position', [180 55 150 20], 'FontSize', 11);
        burstCriterionEdit = uieditfield(badPointPanel, 'numeric', ...
            'Position', [320 55 100 20], ...
            'Value', 20, ...
            'Tag', 'burstCriterionEdit');

        % Window Criterion
        uilabel(badPointPanel, 'Text', 'Window Criterion:', ...
            'Position', [180 30 150 20], 'FontSize', 11);
        windowCriterionEdit = uieditfield(badPointPanel, 'numeric', ...
            'Position', [320 30 100 20], ...
            'Value', 0.25, ...
            'Tag', 'windowCriterionEdit');

        % Window Criterion Tolerances
        uilabel(badPointPanel, 'Text', 'Window Criterion Tolerances:', ...
            'Position', [180 5 200 20], 'FontSize', 11);
        windowToleranceEdit = uieditfield(badPointPanel, 'text', ...
            'Position', [320 5 100 20], ...
            'Value', '[-Inf 7]', ...
            'Tag', 'windowToleranceEdit');

        % Additional instructions
        % uilabel(badPointPanel, 'Text', 'Configure bad point detection criteria.', ...
        %     'Position', [10 250 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);

        % 可以根据需要添加更多参数设置
    end

    function createEpochRestPanel(parent)
        % 创建 "Epoching Resting State Data Parameters" 面板
        epochRestPanel = uipanel(parent, 'Title', 'Epoching Resting State Data Parameters', ...
            'Position', [10 90 620 120], ... % 调整位置以适应滚动
            'BackgroundColor', [0.8 0.9 1], 'FontSize', 12);

        % Epoch Length
        uilabel(epochRestPanel, 'Text', 'Epoch Length (s):', ...
            'Position', [10 60 150 22], 'FontSize', 11);
        epochLengthEdit = uieditfield(epochRestPanel, 'numeric', ...
            'Position', [170 60 100 22], ...
            'Value', 2, ...
            'Tag', 'epochLengthEdit');

        % Epoch Overlap
        uilabel(epochRestPanel, 'Text', 'Epoch Overlap:', ...
            'Position', [10 30 150 22], 'FontSize', 11);
        epochOverlapEdit = uieditfield(epochRestPanel, 'numeric', ...
            'Position', [170 30 100 22], ...
            'Value', 0.5, ...
            'Tag', 'epochOverlapEdit');

        % % Additional instructions
        % uilabel(epochRestPanel, 'Text', 'Configure epoching parameters for resting state data.', ...
        %     'Position', [10 160 600 22], 'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
    end

    %% 生成 params 结构的回调函数
    function generateParamsCallback(fig)
        % 初始化参数结构
        params = struct();

        %% 收集 Channel Information
        chanPanel = findobj(fig, 'Title', 'Channel Information');
        allChanStr = findobj(chanPanel, 'Tag', 'chanAllEdit').Value;
        try
            params.ChanInfo.AllChan = eval(allChanStr);
        catch
            uialert(fig, 'Invalid format for All Channels. Please enter a valid vector (e.g., [1:129]).', 'Input Error');
            return;
        end
        params.ChanInfo.RefChan = findobj(chanPanel, 'Tag', 'refChanEdit').Value;
        eogIdxStr = findobj(chanPanel, 'Tag', 'eogIdxEdit').Value;
        try
            params.ChanInfo.EOGChanIdx = eval(eogIdxStr);
        catch
            uialert(fig, 'Invalid format for EOG Channel Index. Please enter a valid vector (e.g., [1 2 3]).', 'Input Error');
            return;
        end
        params.ChanInfo.EOGChanLabels = findobj(chanPanel, 'Tag', 'eogLabelsEdit').Value;
        knownBadIdxStr = findobj(chanPanel, 'Tag', 'knownBadIdxEdit').Value;
        try
            params.ChanInfo.KnownBadChanIdx = eval(knownBadIdxStr);
        catch
            uialert(fig, 'Invalid format for Known Bad Channel Index. Please enter a valid vector (e.g., [1 2 3]).', 'Input Error');
            return;
        end
        params.ChanInfo.KnownBadChanLabels = findobj(chanPanel, 'Tag', 'knownBadLabelsEdit').Value;

        %% 收集 Downsampling and Filtering Parameters
        downPanel = findobj(fig, 'Title', 'Downsampling and Filtering Parameters');
        params.DownsamplingRate = findobj(downPanel, 'Tag', 'downRateEdit').Value;
        params.Filter.LowCutoff = findobj(downPanel, 'Tag', 'lowCutEdit').Value;
        params.Filter.HighCutoff = findobj(downPanel, 'Tag', 'highCutEdit').Value;
        params.Filter.PowerLineFrequency = findobj(downPanel, 'Tag', 'plfEdit').Value;

        %% 收集 Epoching Parameters
        epochPanel = findobj(fig, 'Title', 'Epoching Parameters');
        params.EpochTask.Marker = findobj(epochPanel, 'Tag', 'markerEdit').Value;
        epochWindowStr = findobj(epochPanel, 'Tag', 'epochWindowEdit').Value;
        try
            params.EpochTask.EpochWindow = eval(epochWindowStr);
        catch
            uialert(fig, 'Invalid format for Epoch Window. Please enter a valid vector (e.g., [1 39]).', 'Input Error');
            return;
        end

        %% 收集 Bad Channel Detection Parameters
        badChanPanel = findobj(fig, 'Title', 'Bad Channel Detection Parameters');
        % Reject Channel Normalization
        rejectNormChk = findobj(badChanPanel, 'Tag', 'rejectNormChk').Value;
        if rejectNormChk
            params.BadChan.rejectchan_NormOn = 'on';
        else
            params.BadChan.rejectchan_NormOn = 'off';
        end

        % Specific Frequency Range
        specFreqChk = findobj(badChanPanel, 'Tag', 'specFreqChk').Value;
        if specFreqChk
            specFreqStr = findobj(badChanPanel, 'Tag', 'specFreqEdit').Value;
            try
                params.BadChan.rejectchan_SpecFreqRange = eval(specFreqStr);
                params.BadChan.rejectchan_SpecOn = 'on';
            catch
                uialert(fig, 'Invalid format for Specific Frequency Range. Please enter a valid vector (e.g., [0.1 40]).', 'Input Error');
                return;
            end
        else
            params.BadChan.rejectchan_SpecFreqRange = [];
            params.BadChan.rejectchan_SpecOn = 'off';
        end

        % Specific Threshold
        specThreshChk = findobj(badChanPanel, 'Tag', 'specThreshChk').Value;
        if specThreshChk
            specThreshStr = findobj(badChanPanel, 'Tag', 'specThreshEdit').Value;
            try
                params.BadChan.rejectchan_SpecThreshold = eval(specThreshStr);
            catch
                uialert(fig, 'Invalid format for Specific Threshold. Please enter a valid vector (e.g., [-5 5]).', 'Input Error');
                return;
            end
        else
            params.BadChan.rejectchan_SpecThreshold = [];
        end

        % Reject Channel Threshold
        rejectThreshStr = findobj(badChanPanel, 'Tag', 'rejectThreshEdit').Value;
        try
            params.BadChan.rejectchan_Threshold = eval(rejectThreshStr);
        catch
            uialert(fig, 'Invalid format for Reject Channel Threshold. Please enter a valid vector (e.g., [-5 5]).', 'Input Error');
            return;
        end

        % Kurtosis Based Rejection
        kurtOnChk = findobj(badChanPanel, 'Tag', 'kurtOnChk').Value;
        if kurtOnChk
            params.BadChan.rejectchan_KurtOn = 'on';
            kurtThreshStr = findobj(badChanPanel, 'Tag', 'kurtThreshEdit').Value;
            try
                params.BadChan.rejectchan_KurtThreshold = eval(kurtThreshStr);
            catch
                uialert(fig, 'Invalid format for Kurtosis Threshold. Please enter a valid vector (e.g., [-5 5]).', 'Input Error');
                return;
            end
        else
            params.BadChan.rejectchan_KurtOn = 'off';
            params.BadChan.rejectchan_KurtThreshold = [];
        end

        % Probability Based Rejection
        probOnChk = findobj(badChanPanel, 'Tag', 'probOnChk').Value;
        if probOnChk
            params.BadChan.rejectchan_ProbOn = 'on';
            probThreshStr = findobj(badChanPanel, 'Tag', 'probThreshEdit').Value;
            try
                params.BadChan.rejectchan_ProbThreshold = eval(probThreshStr);
            catch
                uialert(fig, 'Invalid format for Probability Threshold. Please enter a valid vector (e.g., [-5 5]).', 'Input Error');
                return;
            end
        else
            params.BadChan.rejectchan_ProbOn = 'off';
            params.BadChan.rejectchan_ProbThreshold = [];
        end

        % FASTER Channel Variance
        fasterVarChk = findobj(badChanPanel, 'Tag', 'fasterVarChk').Value;
        if fasterVarChk
            params.BadChan.FASTER_ChanVarOn = 'on';
            params.BadChan.FASTER_ChanVarThreshold = findobj(badChanPanel, 'Tag', 'fasterVarThreshEdit').Value;
        else
            params.BadChan.FASTER_ChanVarOn = 'off';
            params.BadChan.FASTER_ChanVarThreshold = [];
        end

        % FASTER Channel Mean Correlation
        fasterMeanCorrChk = findobj(badChanPanel, 'Tag', 'fasterMeanCorrChk').Value;
        if fasterMeanCorrChk
            params.BadChan.FASTER_ChanMeanCorrOn = 'on';
            params.BadChan.FASTER_ChanMeanCorrThreshold = findobj(badChanPanel, 'Tag', 'fasterMeanCorrThreshEdit').Value;
        else
            params.BadChan.FASTER_ChanMeanCorrOn = 'off';
            params.BadChan.FASTER_ChanMeanCorrThreshold = [];
        end

        % FASTER Channel Hurst
        fasterHurstChk = findobj(badChanPanel, 'Tag', 'fasterHurstChk').Value;
        if fasterHurstChk
            params.BadChan.FASTER_ChanHurst = 'on';
            params.BadChan.FASTER_ChanHurstThreshold = findobj(badChanPanel, 'Tag', 'fasterHurstThreshEdit').Value;
        else
            params.BadChan.FASTER_ChanHurst = 'off';
            params.BadChan.FASTER_ChanHurstThreshold = [];
        end

        % Clean Raw Data Parameters
        % Flat Line Detection
        cleanFlatLineChk = findobj(badChanPanel, 'Tag', 'cleanFlatLineChk').Value;
        if cleanFlatLineChk
            params.BadChan.CleanRawData_FlatLineOn = 'on';
            params.BadChan.CleanRawData_FlatLineThreshold = findobj(badChanPanel, 'Tag', 'cleanFlatLineThreshEdit').Value;
        else
            params.BadChan.CleanRawData_FlatLineOn = 'off';
            params.BadChan.CleanRawData_FlatLineThreshold = [];
        end

        % Clean Drifts High Pass
        cleanDriftsHighPassStr = findobj(badChanPanel, 'Tag', 'cleanDriftsHighPassEdit').Value;
        try
            params.BadChan.CleanRawData_CleanDriftsHighPass = eval(cleanDriftsHighPassStr);
        catch
            uialert(fig, 'Invalid format for Clean Drifts High Pass. Please enter a valid vector (e.g., [0.25, 0.75]).', 'Input Error');
            return;
        end

        % Clean Channel
        cleanChanChk = findobj(badChanPanel, 'Tag', 'cleanChanChk').Value;
        if cleanChanChk
            params.BadChan.CleanRawData_CleanChan = 'on';
        else
            params.BadChan.CleanRawData_CleanChan = 'off';
        end

        % Line Noise Threshold
        params.BadChan.CleanRawData_LineNoiseThreshold = findobj(badChanPanel, 'Tag', 'cleanLineNoiseThreshEdit').Value;

        % Correlation Threshold
        params.BadChan.CleanRawData_CorrThreshold = findobj(badChanPanel, 'Tag', 'cleanCorrThreshEdit').Value;

        %% 收集 Re-referencing Parameters
        reRefPanel = findobj(fig, 'Title', 'Re-referencing Parameters');
        % Channel to Average Checkbox
        chanAvgChk = findobj(reRefPanel, 'Tag', 'chanAvgChk').Value;
        if chanAvgChk
            params.Refer.Chan2averageOn = 'on';
            chanAvgStr = findobj(reRefPanel, 'Tag', 'chanAvgEdit').Value;
            try
                params.Reeref.Chan2average = eval(chanAvgStr);
            catch
                uialert(fig, 'Invalid format for Channels to Average. Please enter a valid vector (e.g., [1 2 3]).', 'Input Error');
                return;
            end
        else
            params.Refer.Chan2averageOn = 'off';
            params.Reeref.Chan2average = [];
        end

        %% 收集 ICA Parameters
        icaPanel = findobj(fig, 'Title', 'ICA Parameters');
        % Filter ICA On
        filterICAonChk = findobj(icaPanel, 'Tag', 'filterICAonChk').Value;
        if filterICAonChk
            params.BadIC.filterICAon = 'on';
        else
            params.BadIC.filterICAon = 'off';
        end

        % ICA Low Cutoff
        params.BadIC.filterICAlocutoff = findobj(icaPanel, 'Tag', 'icaLowCutEdit').Value;

        % IC Label On
        icLabelChk = findobj(icaPanel, 'Tag', 'icLabelChk').Value;
        if icLabelChk
            params.BadIC.ICLabelOn = 'on';
        else
            params.BadIC.ICLabelOn = 'off';
        end

        % FASTER On
        fasterChk = findobj(icaPanel, 'Tag', 'fasterChk').Value;
        if fasterChk
            params.BadIC.FASTEROn = 'on';
        else
            params.BadIC.FASTEROn = 'off';
        end

        % IC Label Threshold
        icLabelThreshStr = findobj(icaPanel, 'Tag', 'icLabelThreshEdit').Value;
        try
            params.BadIC.ICLabelThreshold = eval(icLabelThreshStr);
        catch
            uialert(fig, 'Invalid format for IC Label Threshold. Please enter a valid matrix.', 'Input Error');
            return;
        end

        % Number of Runs
        params.BadIC.Nrun = findobj(icaPanel, 'Tag', 'nRunEdit').Value;

        %% 收集 Bad Point Detection Parameters
        badPointPanel = findobj(fig, 'Title', 'Bad Point Detection Parameters');
        % Bad Point Detection On
        badPointOnChk = findobj(badPointPanel, 'Tag', 'badPointOnChk').Value;
        if badPointOnChk
            params.BadPoint.on = 'on';
            params.BadPoint.BurstCriterion = findobj(badPointPanel, 'Tag', 'burstCriterionEdit').Value;
            params.BadPoint.WindowCriterion = findobj(badPointPanel, 'Tag', 'windowCriterionEdit').Value;
            windowTolStr = findobj(badPointPanel, 'Tag', 'windowToleranceEdit').Value;
            try
                params.BadPoint.WindowCriterionTolerances = eval(windowTolStr);
            catch
                uialert(fig, 'Invalid format for Window Criterion Tolerances. Please enter a valid vector (e.g., [-Inf 7]).', 'Input Error');
                return;
            end
        else
            params.BadPoint.on = 'off';
            params.BadPoint.BurstCriterion = [];
            params.BadPoint.WindowCriterion = [];
            params.BadPoint.WindowCriterionTolerances = [];
        end

        %% 收集 Epoching Resting State Data Parameters
        epochRestPanel = findobj(fig, 'Title', 'Epoching Resting State Data Parameters');
        params.EpochRest.EpochLength = findobj(epochRestPanel, 'Tag', 'epochLengthEdit').Value;
        params.EpochRest.EpochOverlap = findobj(epochRestPanel, 'Tag', 'epochOverlapEdit').Value;

        %% 显示生成的 params 结构
        paramsFig = uifigure('Name', 'Generated Params', 'Position', [1400 100 600 900], ...
                            'Color', [0.8 0.9 1]);
        uitextarea(paramsFig, 'Value', {structToString(params)}, ...
            'Editable', 'off', 'Position', [10 10 580 880], ...
            'FontName', 'Courier New', 'FontSize', 10);

        % 可选：将 params 赋值到基础工作区
        assignin('base', 'params', params);
        delete(fig);
    end

    %% 辅助函数：将结构体转换为字符串
    function str = structToString(s, indent)
        if nargin < 2
            indent = '';
        end
        fields = fieldnames(s);
        str = '';
        for i = 1:length(fields)
            field = fields{i};
            value = s.(field);
            if isstruct(value)
                str = [str, sprintf('%s.%s = struct();\n', indent, field)];
                str = [str, structToString(value, [indent '    '])];
            elseif isnumeric(value) || islogical(value)
                str = [str, sprintf('%s.%s = %s;\n', indent, field, mat2str(value))];
            elseif ischar(value) || isstring(value)
                str = [str, sprintf('%s.%s = ''%s'';\n', indent, field, value)];
            else
                str = [str, sprintf('%s.%s = []; %% Unsupported type\n', indent, field)];
            end
        end
    end
end
