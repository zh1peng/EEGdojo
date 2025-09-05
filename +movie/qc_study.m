function [quality_data_table, qc_table] = qc_study(varargin)
% GET_QC Generate and visualize data quality metrics from EEGLAB .set files.
%   This function scans a directory for .set files, calculates a range of data
%   quality metrics for each, and then flags outliers based on a robust
%   statistical method (Median Absolute Deviation - MAD).
% 
% Input Arguments (Name-Value Pairs):
%   'study_path' (char, REQUIRED)
%       The absolute path to the root directory containing .set files.
% 
%   'searchstring' (char, default: '.*\.set$')
%       A regular expression to find specific .set files within the study_path.
% 
%   'subject_parser' (char, default: '(?<sub>.+)')
%       A regex to extract a subject ID from the filename.
% 
%   'subjects_to_include' ('all' or numeric vector, default: 'all')
%       Specify a subset of subjects to process based on their sorted index.
% 
%   'cutoff_mad' (numeric, default: 3)
%       The number of Median Absolute Deviations (MAD) from the median to
%       use as the cutoff for flagging outliers.
% 
%   'metrics_to_analyze' (cellstr, default: all non-NaN metrics)
%       A cell array of strings specifying which quality metrics to analyze.
% 
%   'save_path' (char, default: '')
%       If a full path is provided, the generated figure is saved.
% 
% Output Arguments:
%   qc_table    - A table where rows are subjects and columns are quality
%                 metrics, containing pass/fail flags (1/0) and an overall
%                 `include_auto` flag.
% 
% Example:
%   qc_table = movie.get_qc('study_path', 'C:\EEG_Data\MovieStudy', 'cutoff_mad', 2.5);

%% 1. Parse Inputs
p = inputParser;
p.FunctionName = 'movie.get_qc';
addParameter(p, 'study_path', '', @ischar);
addParameter(p, 'searchstring', '.*\.set$', @ischar);
addParameter(p, 'subject_parser', '(?<sub>.+)', @ischar);
addParameter(p, 'subjects_to_include', 'all', @(x) (ischar(x) && strcmpi(x, 'all')) || (isnumeric(x) && isvector(x)));
addParameter(p, 'cutoff_mad', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'metrics_to_analyze', {}, @iscellstr);
addParameter(p, 'save_path', '', @ischar);
addParameter(p, 'thresholds', struct(), @isstruct);
parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.study_path), 'You must provide a ''study_path''.');

%% 2. Generate Quality Data Table
fprintf('Study path provided. Generating quality data from .set files...\n');
[paths, names] = filesearch_regexp(opt.study_path, opt.searchstring, true);
assert(~isempty(names), 'No .set files matched searchstring in study_path.');

if exist('natsort', 'file'), [~, sort_idx] = natsort(names); paths = paths(sort_idx); names = names(sort_idx); end
setFiles = fullfile(paths, names);

% --- Apply subject selection ---
if isnumeric(opt.subjects_to_include)
    num_found = numel(setFiles);
    indices = opt.subjects_to_include;
    indices = unique(round(indices(indices > 0 & indices <= num_found)));
    assert(~isempty(indices), 'The specified subjects_to_include indices resulted in an empty list of subjects.');
    fprintf('Including %d subjects based on specified indices.\n', numel(indices));
    setFiles = setFiles(indices);
end

quality_reports = {};
for i = 1:numel(setFiles)
    fpath = setFiles{i};
    [subID, ~] = parseSubSesFromFile(fpath, opt.subject_parser);
    subKey = makeFieldKey(subID);
    
    fprintf('  Loading and calculating quality for %s...\n', subID);
    EEG = pop_loadset(fpath);
    if ~isempty(EEG.data)
            % Calculate data quality metrics
            EEG2qc =pop_eegfiltnew(EEG, 'locutoff', 1);
            EEG2qc = eeg_checkset(EEG2qc);
            data_quality = calculate_data_quality(EEG2qc);
    else
            data_quality = struct();
    end
    report_row = struct2table(data_quality, 'AsArray', true);
    report_row.subject = string(subKey);
    quality_reports{end+1} = report_row;
end
quality_data_table = vertcat(quality_reports{:});

%% 3. Determine which metrics to analyze and plot
if ~ismember('subject', quality_data_table.Properties.VariableNames), error('Quality data table must contain a ''subject'' column.'); end
if ~strcmp(quality_data_table.Properties.VariableNames{1}, 'subject'), quality_data_table = movevars(quality_data_table, 'subject', 'Before', 1); end

num_subjects = height(quality_data_table);
assert(num_subjects > 0, 'Quality data table is empty.');
quality_data_struct = table2struct(quality_data_table, 'ToScalar', true);

if isempty(opt.metrics_to_analyze)
    metrics_to_analyze = {};
    all_possible_metrics = quality_data_table.Properties.VariableNames(2:end);
    for i = 1:numel(all_possible_metrics)
        metric_name = all_possible_metrics{i};
        if iscell(quality_data_struct.(metric_name)) || all(isnan(quality_data_struct.(metric_name))), continue; end
        metrics_to_analyze{end+1} = metric_name;
    end
else
    metrics_to_analyze = opt.metrics_to_analyze;
end

if isempty(metrics_to_analyze), warning('No valid quality metrics found to analyze.'); qc_table = table(); return; end

%% 4. Calculate Cutoffs and Generate Flags
qc_table = table();
qc_table.subject = quality_data_table.subject;
cutoffs = struct();

higher_is_worse = {'percent_bad_channels', 'percent_rejected_epochs', 'global_field_power_std', ...
                   'amplitude_spread_proxy', 'high_freq_noise_ratio', 'line_noise_power_50Hz', ...
                   'line_noise_power_60Hz', 'percent_rejected_by_windows'};

for i = 1:numel(metrics_to_analyze)
    metric = metrics_to_analyze{i};
    vals = quality_data_struct.(metric);
    
    if isfield(opt.thresholds, metric)
        cutoff = opt.thresholds.(metric);
        fprintf('Using user-defined threshold for %s: %.2f\n', metric, cutoff);
    else
        med = custom_nanmedian(vals);
        mad_val = custom_nanmad(vals);
        
        if any(strcmpi(metric, higher_is_worse))
            cutoff = med + opt.cutoff_mad * mad_val;
        else
            cutoff = med - opt.cutoff_mad * mad_val;
        end
    end

    if any(strcmpi(metric, higher_is_worse))
        pass_fail_flags = vals <= cutoff;
    else
        pass_fail_flags = vals >= cutoff;
    end
    
    pass_fail_flags(isnan(vals)) = true;
    qc_table.(metric) = pass_fail_flags;
    cutoffs.(metric) = cutoff;
end

metric_flags = table2array(qc_table(:, 2:end));
qc_table.include_auto = all(metric_flags, 2);

%% 5. Create Visualization
num_plots = numel(metrics_to_analyze);
if num_plots == 0, return; end

figure('Name', 'Data Quality Check', 'Color', 'w', 'Position', [100 100 1400 400*ceil(num_plots/4)]);
sgtitle(sprintf('Data Quality Metrics (Cutoff = median +/- %.1f*MAD)', opt.cutoff_mad), 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:numel(metrics_to_analyze)
    metric = metrics_to_analyze{i};
    vals = quality_data_struct.(metric);
    flags = qc_table.(metric);
    
    if num_plots == 1
        subplot(1, 1, 1);
    else
        subplot(ceil(num_plots/4), 4, i);
    end
    
    pass_idx = find(flags == 1);
    fail_idx = find(flags == 0);
    
    hold on;
    if ~isempty(pass_idx), bar(pass_idx, vals(pass_idx), 'FaceColor', [0 .6 0.2], 'EdgeColor', 'none'); end
    if ~isempty(fail_idx), bar(fail_idx, vals(fail_idx), 'FaceColor', [0.8 0.2 0.1], 'EdgeColor', 'none'); end
    
    line_h = yline(cutoffs.(metric), '--r', 'LineWidth', 1.5);

    current_ylim = ylim;
    new_ylim_top = current_ylim(2) * 1.2;
    ylim([current_ylim(1), new_ylim_top]);

    legend(line_h, sprintf('Cutoff (%.2f)', cutoffs.(metric)), 'Location', 'northeast', 'Box', 'off');
    hold off;
    
    title(strrep(metric, '_', ' '), 'FontWeight', 'normal');
    ylabel('Value');
    xlabel('Subject');
    xlim([0 num_subjects+1]);
    % set(gca, 'XTick', 1:num_subjects, 'XTickLabel', strrep(quality_data_table.subject, 'sub_', ''), 'XTickLabelRotation', 90);
    % set Xtic off
    set(gca, 'XTickLabel', []);
    grid on;
end

%% 6. Save figure if requested
if ~isempty(opt.save_path)
    fprintf('Saving quality report figure to: %s\n', opt.save_path);
    exportgraphics(gcf, opt.save_path, 'Resolution', 150);
end

end

%% ----------------- Local helpers -----------------
function [subID, sesID] = parseSubSesFromFile(fpath, parserPattern)
    [~,fname,~] = fileparts(fpath);
    m = regexp(fname, parserPattern, 'names');
    assert(~isempty(m) && isfield(m,'sub'), 'Failed to parse subject from: %s', fname);
    subID = m.sub;
    if isfield(m,'ses') && ~isempty(m.ses), sesID = m.ses; else, sesID = ''''; end
end

function key = makeFieldKey(txt)
    key = lower(regexprep(txt,'[^A-Za-z0-9]','_'));
    if ~isempty(key) && isstrprop(key(1),'digit'), key = ['x' key]; end
end

function m = custom_nanmedian(x)
    x = x(~isnan(x));
    if isempty(x)
        m = NaN;
        return;
    end
    x = sort(x);
    n = length(x);
    if mod(n, 2) == 1
        m = x((n + 1) / 2);
    else
        m = (x(n / 2) + x(n / 2 + 1)) / 2;
    end
end

function m = custom_nanmad(x)
    med = custom_nanmedian(x);
    m = custom_nanmedian(abs(x - med));
end
