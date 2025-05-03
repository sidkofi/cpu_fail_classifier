% VISUALIZATIONS
%
% Description:
%   This script generates all the MATLAB figures used in the accompanying
%   README.md file. It should be run AFTER the main
%   classification script (e.g., classify_cpu_failures.m) has been executed,
%   as it relies on variables created by that script being present in the
%   MATLAB workspace.
%
% Prerequisites (must exist in workspace):
%   - data: Table loaded from benchmark_data.csv
%   - new_features: Matrix of calculated features [pk_count, mag_ratio, pk_dist_sep, mad_val]
%   - true_labels_valid: Cell array of true labels for valid runs
%   - predicted_labels_valid: Cell array of predicted labels for valid runs
%   - valid_indices: Logical array indicating valid runs
%   - Fs, ma_window, butter_order, butter_cutoff, wavelet_name, num_levels
%   - peak_height_factor_sig, peak_dist_sec_sig, peak_exclude_window_samples
%   - min_dist_for_feature_sec, min_peak_height_factor_for_dist
%   - threshold_mad_max_for_none, threshold_pk_count_min_for_none,
%     threshold_freq_osc_min_peak_count, threshold_freq_osc_max_ratio,
%     threshold_mad_min_for_freq_osc, threshold_thermal_min_ratio_final,
%     threshold_mad_min_for_thermal
%
% Output:
%   - Saves multiple PNG files corresponding to the figures needed for the
%     blog post into the 'images' subdirectory (relative to current folder).
%     Ensure the 'generated_images' folder exists before running.

% --- Create Output Directory ---
output_dir = 'generated_images';
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
   fprintf('Created output directory: %s\n', output_dir);
end

% --- Check for necessary variables ---
required_vars = {'data', 'new_features', 'true_labels_valid', 'predicted_labels_valid', ...
                 'valid_indices', 'Fs', 'ma_window', 'butter_order', 'butter_cutoff', ...
                 'wavelet_name', 'num_levels', 'peak_height_factor_sig', 'peak_dist_sec_sig', ...
                 'peak_exclude_window_samples', 'min_dist_for_feature_sec', ...
                 'min_peak_height_factor_for_dist', 'threshold_mad_max_for_none', ...
                 'threshold_pk_count_min_for_none', 'threshold_freq_osc_min_peak_count', ...
                 'threshold_freq_osc_max_ratio', 'threshold_mad_min_for_freq_osc', ...
                 'threshold_thermal_min_ratio_final', 'threshold_mad_min_for_thermal'};

missing_vars = {};
for k=1:length(required_vars)
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', required_vars{k}))
        missing_vars{end+1} = required_vars{k};
    end
end
if ~isempty(missing_vars)
    error('Missing required variables in workspace: %s. Run main classification script first.', strjoin(missing_vars, ', '));
end
fprintf('All required variables found in workspace.\n');

% --- Define Example Runs ---
fprintf('\n--- Defining Example Runs for Visualization ---\n');
example_runs = struct();
try
    failure_types_unique = categories(unique(categorical(true_labels_valid)));
    original_run_indices_map = find(valid_indices);

    for i = 1:length(failure_types_unique)
        type = failure_types_unique{i};
        idx_in_valid = find(strcmp(true_labels_valid, type), 1);
        if ~isempty(idx_in_valid)
            example_runs.(type) = original_run_indices_map(idx_in_valid);
            fprintf('Selected Run ID %d as example for %s\n', example_runs.(type), type);
        else
            warning('Could not find a valid example run for type: %s', type);
            example_runs.(type) = -1;
        end
    end
catch ME_example
    error('Could not automatically define example runs. Error: %s', ME_example.message);
end

%% Plot 1: Example Raw Signals (Separate Figures, Dynamic Y-Limits)
fprintf('\n--- Plotting Example Raw Signals ---\n');
failure_types_to_plot_raw = fieldnames(example_runs)';
num_types_raw = length(failure_types_to_plot_raw);

for i = 1:num_types_raw
    current_type_str = failure_types_to_plot_raw{i};
    example_run_id = example_runs.(current_type_str);
    if example_run_id == -1, continue; end

    run_subset = data(data.Run_ID == example_run_id, :);
    if isempty(run_subset) || ~ismember('Time', run_subset.Properties.VariableNames) || ~ismember('Frequency_GHz', run_subset.Properties.VariableNames)
        warning('Example run %d (%s) invalid for plotting raw signal.', example_run_id, current_type_str); continue;
    end
    time_data = run_subset.Time;
    freq_ghz_data = run_subset.Frequency_GHz;
    if isempty(time_data) || isempty(freq_ghz_data) || length(time_data) ~= length(freq_ghz_data) || all(isnan(freq_ghz_data))
        warning('Invalid data for Run ID %d (%s) raw plot.', example_run_id, current_type_str); continue;
    end

    min_freq_run = min(freq_ghz_data); max_freq_run = max(freq_ghz_data);
    y_padding = (max_freq_run - min_freq_run) * 0.05; if y_padding == 0, y_padding = 0.1; end
    ylim_dynamic = [min_freq_run - y_padding, max_freq_run + y_padding];

    fig_raw = figure('Name', sprintf('Example Raw: %s', current_type_str), 'Position', [100 + (i-1)*40, 100 + (i-1)*40, 600, 400], 'Visible', 'off');
    plot(time_data, freq_ghz_data, 'LineWidth', 1.5); grid on; box on;
    ylim(ylim_dynamic); title(strrep(current_type_str,'_',' '), 'FontSize', 12);
    xlabel('Time (s)', 'FontSize', 10); ylabel('Frequency (GHz)', 'FontSize', 10); set(gca, 'FontSize', 10);
    filename_raw = fullfile(output_dir, sprintf('fail_mode_%s.png', lower(current_type_str)));
    try
        saveas(fig_raw, filename_raw);
        fprintf('Saved: %s\n', filename_raw);
    catch ME_save
        warning('Could not save figure %s. Error: %s', filename_raw, ME_save.message);
    end
    close(fig_raw);
end

%% Plot 2: Filter Responses
fprintf('\n--- Plotting Filter Responses ---\n');
% Moving Average
fig_ma = figure('Name','Filter Response: Moving Average', 'Visible', 'off');
ma_coeffs = ones(1, ma_window) / ma_window; [h_ma, w_ma] = freqz(ma_coeffs, 1, 1024); f_ma = w_ma/(2*pi) * Fs;
plot(f_ma, 20*log10(abs(h_ma)), 'LineWidth', 2); title('Frequency Response'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on; box on; xlim([0, Fs/2]); ylim([-30, 5]);
ma_params_str = sprintf('Window Size: %d', ma_window);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ma_params_str, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FitBoxToText', 'on');
filename_ma = fullfile(output_dir, 'filter_response_ma.png');
try saveas(fig_ma, filename_ma); fprintf('Saved: %s\n', filename_ma); catch ME_save, warning('Could not save figure %s. Error: %s', filename_ma, ME_save.message); end
close(fig_ma);

% Butterworth
fig_butter = figure('Name','Filter Response: Butterworth', 'Visible', 'off');
[b, a] = butter(butter_order, butter_cutoff/(Fs/2)); [h_butter, w_butter] = freqz(b, a, 1024); f_butter = w_butter/(2*pi) * Fs;
plot(f_butter, 20*log10(abs(h_butter)), 'LineWidth', 2); title('Frequency Response'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on; box on; xlim([0, Fs/2]); ylim([-30, 5]);
butter_params_str = sprintf('Order: %d Cutoff: %.1f Hz', butter_order, butter_cutoff);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', butter_params_str, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FitBoxToText', 'on');
filename_butter = fullfile(output_dir, 'filter_response_butter.png');
try saveas(fig_butter, filename_butter); fprintf('Saved: %s\n', filename_butter); catch ME_save, warning('Could not save figure %s. Error: %s', filename_butter, ME_save.message); end
close(fig_butter);

%% Plot 3: Wavelet Functions (Saving Phi and Psi separately)
fprintf('\n--- Plotting Wavelet Functions ---\n');
iterations = 8;
try
    [phi, psi, xval] = wavefun(wavelet_name, iterations);

    % Plot and Save Scaling Function (phi)
    fig_phi = figure('Name', sprintf('%s Scaling Function', wavelet_name), 'Position', [100, 100, 350, 350], 'Visible', 'off'); % Smaller figure
    plot(xval, phi, 'LineWidth', 1.5); title('Scaling Function (\phi)'); xlabel('Time'); ylabel('Amplitude'); grid on; box on; axis tight;
    filename_phi = fullfile(output_dir, sprintf('wavelet_%s_phi.png', wavelet_name));
    try saveas(fig_phi, filename_phi); fprintf('Saved: %s\n', filename_phi); catch ME_save, warning('Could not save figure %s. Error: %s', filename_phi, ME_save.message); end
    close(fig_phi);

    % Plot and Save Wavelet Function (psi)
    fig_psi = figure('Name', sprintf('%s Wavelet Function', wavelet_name), 'Position', [150, 150, 350, 350], 'Visible', 'off'); % Smaller figure
    plot(xval, psi, 'LineWidth', 1.5); title('Wavelet Function (\psi)'); xlabel('Time'); ylabel('Amplitude'); grid on; box on; axis tight;
    filename_psi = fullfile(output_dir, sprintf('wavelet_%s_psi.png', wavelet_name));
    try saveas(fig_psi, filename_psi); fprintf('Saved: %s\n', filename_psi); catch ME_save, warning('Could not save figure %s. Error: %s', filename_psi, ME_save.message); end
    close(fig_psi);

catch ME_wavefun
    warning('Could not plot wavelet functions. Wavelet Toolbox installed? Error: %s', ME_wavefun.message);
end

%% Plot 4: D1/D2 Reconstructions (Separate figure per type)
fprintf('\n--- Plotting D1 and D2 Reconstructions ---\n');
focus_levels_d1d2 = [1, 2];
failure_types_to_plot_d1d2 = fieldnames(example_runs)';

for i = 1:length(failure_types_to_plot_d1d2)
    type = failure_types_to_plot_d1d2{i};
    idx = example_runs.(type);
    if idx == -1, continue; end

    run_subset = data(data.Run_ID == idx, :);
    if isempty(run_subset) || ~ismember('Time', run_subset.Properties.VariableNames) || ~ismember('Frequency_GHz', run_subset.Properties.VariableNames)
        warning('Run ID %d (%s) invalid for D1/D2 plot.', idx, type); continue; end
    time_data = run_subset.Time; raw_freq_data = run_subset.Frequency_GHz;
    if isempty(time_data) || isempty(raw_freq_data) || length(time_data) ~= length(raw_freq_data) || all(isnan(raw_freq_data))
        warning('Invalid data for Run ID %d (%s) D1/D2 plot.', idx, type); continue; end

    freq_ma = movmean(raw_freq_data, ma_window, 'omitnan'); [b_filt, a_filt] = butter(butter_order, butter_cutoff/(Fs/2)); freq_data = filtfilt(b_filt, a_filt, freq_ma);
    try
        [c, l] = wavedec(freq_data, num_levels, wavelet_name);
    catch ME_wavedec
         warning('Wavedec failed for Run ID %d (%s). Skipping D1/D2 plot. Error: %s', idx, type, ME_wavedec.message); continue;
    end

    fig_d1d2 = figure('Name', sprintf('%s - D1 vs D2', strrep(type, '_', ' ')), 'Position', [100, 100, 900, 500], 'Visible', 'off');

    for i_plot = 1:length(focus_levels_d1d2)
        level = focus_levels_d1d2(i_plot);
        subplot(2, 1, i_plot);
        try
            detail_recon = wrcoef('d', c, l, wavelet_name, level);
            plot(time_data, detail_recon, 'LineWidth', 1.5);
            min_peak_h = std(abs(detail_recon)) * 1.5; [pks, locs] = findpeaks(abs(detail_recon), 'MinPeakHeight', min_peak_h); num_peaks = length(pks);
            if ~isempty(pks) && ~isempty(locs)
                valid_locs = locs(locs <= length(time_data));
                if ~isempty(valid_locs), hold on; plot(time_data(valid_locs), detail_recon(valid_locs), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); hold off; end
            end
            freq_low = (Fs / (2^(level+1))); freq_high = (Fs / (2^level));
            title(sprintf('Detail Level %d Reconstruction (Approx. %.2f-%.2f Hz)', level, freq_low, freq_high)); ylabel('Amplitude'); grid on; box on; xlim([time_data(1), time_data(end)]);
            peak_str = sprintf('Peaks: %d', num_peaks); text(0.98, 0.95, peak_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9, 'BackgroundColor', 'w','EdgeColor','k','Margin', 1);
            if i_plot == length(focus_levels_d1d2), xlabel('Time (s)'); end
        catch ME_detail
             warning('Error plotting detail level %d for Run ID %d (%s). Error: %s', level, idx, type, ME_detail.message);
             title(sprintf('Detail Level %d - Error', level)); % Indicate error on plot
        end
    end
    sgtitle(sprintf('Wavelet Detail Levels for %s (Run %d)', strrep(type, '_', ' '), idx), 'FontSize', 12, 'FontWeight', 'bold');

    filename_d1d2 = fullfile(output_dir, sprintf('wavelet_details_%s.png', lower(type)));
    try saveas(fig_d1d2, filename_d1d2); fprintf('Saved: %s\n', filename_d1d2); catch ME_save, warning('Could not save figure %s. Error: %s', filename_d1d2, ME_save.message); end
    close(fig_d1d2);
end

%% Plot 5: Individual Feature Visualizations
fprintf('\n--- Plotting Individual Feature Visualizations ---\n');

% --- 5a: Peak Count ---
try
    example_type_pk = 'Frequency_Oscillation';
    example_run_id_pk = example_runs.(example_type_pk);
    if example_run_id_pk ~= -1
        run_subset = data(data.Run_ID == example_run_id_pk, :); time_data = run_subset.Time; raw_freq_data = run_subset.Frequency_GHz;
        freq_ma = movmean(raw_freq_data, ma_window, 'omitnan'); [b_filt, a_filt] = butter(butter_order, butter_cutoff/(Fs/2)); freq_data = filtfilt(b_filt, a_filt, freq_ma);
        [c, l] = wavedec(freq_data, num_levels, wavelet_name); detail_recon = wrcoef('d', c, l, wavelet_name, 1); abs_detail_recon = abs(detail_recon);
        min_peak_h = peak_height_factor_sig * std(abs_detail_recon); min_peak_d_samples = round(peak_dist_sec_sig * Fs); if min_peak_d_samples < 1, min_peak_d_samples = 1; end
        [pks, locs] = findpeaks(abs_detail_recon, 'MinPeakHeight', min_peak_h, 'MinPeakDistance', min_peak_d_samples); pk_count_value = length(pks);
        fig_pk = figure('Name', 'Feature Visualization: Peak Count', 'Visible', 'off');
        plot(time_data, abs_detail_recon, 'b-', 'LineWidth', 1, 'DisplayName', '|D1 Signal|'); hold on;
        hline = refline(0, min_peak_h); hline.Color = [0.8 0.2 0.2]; hline.LineStyle = '--'; hline.LineWidth = 1; hline.DisplayName = sprintf('Min Peak Height (%.4f)', min_peak_h);
        if ~isempty(pks) && ~isempty(locs), valid_locs = locs(locs <= length(time_data)); if ~isempty(valid_locs), plot(time_data(valid_locs), pks, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Significant Peaks'); end; end
        hold off; grid on; box on; title(sprintf('Significant Peak Count (Run %d - %s)', example_run_id_pk, strrep(example_type_pk,'_',' '))); xlabel('Time (s)'); ylabel('Absolute D1 Amplitude'); legend('show', 'Location', 'northeast'); xlim([time_data(1), time_data(end)]);
        annotation_str = sprintf('Significant Peak Count: %d', pk_count_value); text(0.98, 0.95, annotation_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 1);
        filename_pk = fullfile(output_dir, 'feature_pk_count_vis.png');
        try saveas(fig_pk, filename_pk); fprintf('Saved: %s\n', filename_pk); catch ME_save, warning('Could not save figure %s. Error: %s', filename_pk, ME_save.message); end
        close(fig_pk);
    else
        warning('Skipping Peak Count visualization: No example run for %s.', example_type_pk);
    end
catch ME_vis, warning('Error generating Peak Count visualization: %s', ME_vis.message); end

% --- 5b: Magnitude Ratio (Separated Peaks) ---
try
    example_type_ratio = 'Thermal_Throttling';
    example_run_id_ratio = example_runs.(example_type_ratio);
     if example_run_id_ratio ~= -1
        run_subset = data(data.Run_ID == example_run_id_ratio, :); time_data = run_subset.Time; raw_freq_data = run_subset.Frequency_GHz; n_samples = length(time_data);
        freq_ma = movmean(raw_freq_data, ma_window, 'omitnan'); [b_filt, a_filt] = butter(butter_order, butter_cutoff/(Fs/2)); freq_data = filtfilt(b_filt, a_filt, freq_ma);
        [c, l] = wavedec(freq_data, num_levels, wavelet_name); detail_recon = wrcoef('d', c, l, wavelet_name, 1); abs_detail_recon = abs(detail_recon);
        mag_ratio_value = 0; avg_top_peak_mag = 0; avg_rest_mag = 0; pks_separated = []; locs_separated = [];
        min_peak_h_ratio = min_peak_height_factor_for_dist * std(abs_detail_recon);
        min_dist_samples_ratio = round(min_dist_for_feature_sec * Fs); if min_dist_samples_ratio < 1, min_dist_samples_ratio = 1; end
        [pks_all, locs_all] = findpeaks(abs_detail_recon, 'MinPeakHeight', min_peak_h_ratio, 'SortStr', 'descend');
        if ~isempty(pks_all)
            pks_separated(1) = pks_all(1); locs_separated(1) = locs_all(1);
            if length(pks_all) >= 2
                first_peak_loc = locs_all(1);
                for j = 2:length(pks_all)
                    second_peak_loc_candidate = locs_all(j);
                    if abs(second_peak_loc_candidate - first_peak_loc) >= min_dist_samples_ratio, pks_separated(2) = pks_all(j); locs_separated(2) = second_peak_loc_candidate; break; end
                end
            end
            avg_top_peak_mag = mean(pks_separated); mask = true(n_samples, 1);
            for p = 1:length(locs_separated)
                 current_loc = locs_separated(p);
                 if current_loc > 0 && current_loc <= n_samples, win_start = max(1, current_loc - peak_exclude_window_samples); win_end = min(n_samples, current_loc + peak_exclude_window_samples); mask(win_start:win_end) = false; end
            end
            if any(mask), avg_rest_mag = mean(abs_detail_recon(mask)); if avg_rest_mag > 1e-9, mag_ratio_value = avg_top_peak_mag / avg_rest_mag; else, mag_ratio_value = avg_top_peak_mag / (1e-9); end
            else, warning('Mask excluded all samples for Run ID %d ratio plot.', example_run_id_ratio); mag_ratio_value = Inf; end
        end
        fig_ratio = figure('Name', 'Feature Visualization: Magnitude Ratio (Separated Peaks)', 'Visible', 'off');
        plot(time_data, abs_detail_recon, 'b-', 'LineWidth', 1, 'DisplayName', '|D1 Signal|'); hold on;
        if avg_rest_mag > 0, hline_rest = refline(0, avg_rest_mag); hline_rest.Color = [0.2 0.8 0.2]; hline_rest.LineStyle = ':'; hline_rest.LineWidth = 1.5; hline_rest.DisplayName = sprintf('Avg Rest Mag (%.4f)', avg_rest_mag); end
        if ~isempty(pks_separated) && ~isempty(locs_separated)
            valid_locs = locs_separated(locs_separated <= length(time_data)); valid_pks = pks_separated(locs_separated <= length(time_data));
            if ~isempty(valid_locs)
                plot(time_data(valid_locs(1)), valid_pks(1), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Peak 1 (Largest)');
                if length(valid_locs) >= 2, plot(time_data(valid_locs(2)), valid_pks(2), 'g*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Peak 2 (Separated)'); end
            end
        end
        hold off; grid on; box on; title(sprintf('Magnitude Ratio (Separated Peaks) (Run %d - %s)', example_run_id_ratio, strrep(example_type_ratio,'_',' '))); xlabel('Time (s)'); ylabel('Absolute D1 Amplitude'); legend('show', 'Location', 'northeast'); xlim([time_data(1), time_data(end)]);
        annotation_str = sprintf('Separated Peak Ratio: %.1f\n(Avg Sep. Peak(s) / Avg Rest)', mag_ratio_value); text(0.98, 0.95, annotation_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 1);
        filename_ratio = fullfile(output_dir, 'feature_mag_ratio_vis.png');
        try saveas(fig_ratio, filename_ratio); fprintf('Saved: %s\n', filename_ratio); catch ME_save, warning('Could not save figure %s. Error: %s', filename_ratio, ME_save.message); end
        close(fig_ratio);
    else
        warning('Skipping Magnitude Ratio visualization: No example run for %s.', example_type_ratio);
    end
catch ME_vis, warning('Error generating Magnitude Ratio visualization: %s', ME_vis.message); end

% --- 5c: Separated Peak Distance ---
try
    example_type_dist = 'Stuck_Frequency';
    example_run_id_dist = example_runs.(example_type_dist);
     if example_run_id_dist ~= -1
        run_subset = data(data.Run_ID == example_run_id_dist, :); time_data = run_subset.Time; raw_freq_data = run_subset.Frequency_GHz;
        freq_ma = movmean(raw_freq_data, ma_window, 'omitnan'); [b_filt, a_filt] = butter(butter_order, butter_cutoff/(Fs/2)); freq_data = filtfilt(b_filt, a_filt, freq_ma);
        [c, l] = wavedec(freq_data, num_levels, wavelet_name); detail_recon = wrcoef('d', c, l, wavelet_name, 1); abs_detail_recon = abs(detail_recon);
        pk_dist_sep_value = 0; min_peak_h_dist = min_peak_height_factor_for_dist * std(abs_detail_recon); min_dist_samples = round(min_dist_for_feature_sec * Fs); if min_dist_samples < 1, min_dist_samples = 1; end
        [pks_all, locs_all] = findpeaks(abs_detail_recon, 'MinPeakHeight', min_peak_h_dist, 'SortStr', 'descend'); locs_separated_peaks = [];
        if length(pks_all) >= 2
            first_peak_loc = locs_all(1); locs_separated_peaks(1) = first_peak_loc;
            for j = 2:length(pks_all)
                second_peak_loc_candidate = locs_all(j);
                if abs(second_peak_loc_candidate - first_peak_loc) >= min_dist_samples, pk_dist_sep_value = abs(second_peak_loc_candidate - first_peak_loc) / Fs; locs_separated_peaks(2) = second_peak_loc_candidate; break; end
            end
        end
        fig_dist = figure('Name', 'Feature Visualization: Separated Peak Distance', 'Visible', 'off');
        plot(time_data, abs_detail_recon, 'b-', 'LineWidth', 1, 'DisplayName', '|D1 Signal|'); hold on;
        plot_handles_dist = []; legend_entries_dist = {};
        if ~isempty(locs_separated_peaks)
            valid_locs_sep = locs_separated_peaks(locs_separated_peaks <= length(time_data));
            if ~isempty(valid_locs_sep)
                pks_to_plot = abs_detail_recon(valid_locs_sep);
                h1 = plot(time_data(valid_locs_sep(1)), pks_to_plot(1), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Peak 1'); plot_handles_dist(end+1) = h1; legend_entries_dist{end+1} = 'Peak 1';
                if length(valid_locs_sep) >= 2
                    h2 = plot(time_data(valid_locs_sep(2)), pks_to_plot(2), 'g*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Peak 2'); plot_handles_dist(end+1) = h2; legend_entries_dist{end+1} = 'Peak 2';
                    y_line = max(pks_to_plot) * 1.1; if y_line == 0, y_line = 0.1*max(abs_detail_recon); end
                    plot([time_data(valid_locs_sep(1)), time_data(valid_locs_sep(2))], [y_line, y_line], 'k--<', 'LineWidth', 1, 'HandleVisibility', 'off'); plot([time_data(valid_locs_sep(1)), time_data(valid_locs_sep(2))], [y_line, y_line], 'k-->', 'LineWidth', 1, 'HandleVisibility', 'off');
                    text(mean(time_data(valid_locs_sep)), y_line, sprintf(' %.1f s ', pk_dist_sep_value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'BackgroundColor', 'w');
                end
            end
        end
        hold off; grid on; box on; title(sprintf('Separated Peak Distance - %s',strrep(example_type_dist,'_',' '))); xlabel('Time (s)'); ylabel('Absolute D1 Amplitude');
        if ~isempty(plot_handles_dist), legend(plot_handles_dist, legend_entries_dist, 'Location', 'best'); else, legend('|D1 Signal|', 'Location', 'best'); end
        xlim([time_data(1), time_data(end)]);
        annotation_str = sprintf('Separated Distance: %.1f s', pk_dist_sep_value); text(0.98, 0.95, annotation_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 1);
        filename_dist = fullfile(output_dir, 'feature_sep_dist_vis.png');
        try saveas(fig_dist, filename_dist); fprintf('Saved: %s\n', filename_dist); catch ME_save, warning('Could not save figure %s. Error: %s', filename_dist, ME_save.message); end
        close(fig_dist);
    else
        warning('Skipping Separated Distance visualization: No example run for %s.', example_type_dist);
    end
catch ME_vis, warning('Error generating Separated Peak Distance visualization: %s', ME_vis.message); end

% --- 5d: Mean Absolute Deviation (MAD) ---
try
    example_type_mad = 'None';
    example_run_id_mad = example_runs.(example_type_mad);
    if example_run_id_mad ~= -1
        run_subset = data(data.Run_ID == example_run_id_mad, :); time_data = run_subset.Time; raw_freq_data = run_subset.Frequency_GHz;
        freq_ma = movmean(raw_freq_data, ma_window, 'omitnan'); [b_filt, a_filt] = butter(butter_order, butter_cutoff/(Fs/2)); freq_data = filtfilt(b_filt, a_filt, freq_ma);
        [c, l] = wavedec(freq_data, num_levels, wavelet_name); detail_recon = wrcoef('d', c, l, wavelet_name, 1); abs_detail_recon = abs(detail_recon);
        mad_value = 0; mean_abs_d1 = 0;
        if ~isempty(abs_detail_recon), mean_abs_d1 = mean(abs_detail_recon); mad_value = mean(abs(abs_detail_recon - mean_abs_d1)); end
        fig_mad = figure('Name', 'Feature Visualization: Mean Absolute Deviation', 'Visible', 'off');
        plot(time_data, abs_detail_recon, 'b-', 'LineWidth', 1, 'DisplayName', '|D1 Signal|'); hold on;
        if mean_abs_d1 >= 0, hline_mean = refline(0, mean_abs_d1); hline_mean.Color = [1 0.5 0]; hline_mean.LineStyle = '--'; hline_mean.LineWidth = 1.5; hline_mean.DisplayName = sprintf('Mean |D1| (%.6f)', mean_abs_d1); end
        hold off; grid on; box on; title(sprintf('Mean Absolute Deviation (MAD) (Run %d - %s)', example_run_id_mad, strrep(example_type_mad,'_',' '))); xlabel('Time (s)'); ylabel('Absolute D1 Amplitude'); legend('show', 'Location', 'northeast'); xlim([time_data(1), time_data(end)]);
        current_ylim = ylim; ylim([min(current_ylim(1), -0.0005) max(current_ylim(2), mean_abs_d1*1.5 + 0.001)]);
        annotation_str = sprintf('MAD: %.6f', mad_value); text(0.98, 0.95, annotation_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 1);
        filename_mad = fullfile(output_dir, 'feature_mad_vis.png');
        try saveas(fig_mad, filename_mad); fprintf('Saved: %s\n', filename_mad); catch ME_save, warning('Could not save figure %s. Error: %s', filename_mad, ME_save.message); end
        close(fig_mad);
    else
        warning('Skipping MAD visualization: No example run for %s.', example_type_mad);
    end
catch ME_vis, warning('Error generating MAD visualization: %s', ME_vis.message); end

%% Plot 6: Confusion Matrix
fprintf('\n--- Plotting Confusion Matrix ---\n');
try
    class_order = {'None', 'Frequency_Oscillation', 'Stuck_Frequency', 'Thermal_Throttling'};
    true_cat = categorical(true_labels_valid, class_order);
    pred_cat = categorical(predicted_labels_valid, class_order);
    fig_cm = figure('Name', 'Confusion Matrix', 'Position', [100 100 550 450], 'Visible', 'off');
    cm = confusionchart(true_cat, pred_cat);
    overall_accuracy = sum(diag(cm.NormalizedValues)) / sum(cm.NormalizedValues(:)) * 100;
    cm.Title = sprintf('CPU Failure Classification (Overall Acc: %.2f%%)', overall_accuracy);
    cm.ColumnSummary = 'column-normalized'; cm.RowSummary = 'row-normalized';
    filename_cm = fullfile(output_dir, 'confusion_matrix.png');
    try saveas(fig_cm, filename_cm); fprintf('Saved: %s\n', filename_cm); catch ME_save, warning('Could not save figure %s. Error: %s', filename_cm, ME_save.message); end
    close(fig_cm);
catch ME_cm
    warning('Error generating Confusion Matrix plot: %s', ME_cm.message);
end
fprintf('\n--- All Visualizations Generated ---\n');
