function generate_cpu_benchmark_data(num_runs)
% Creates synthetic CPU benchmark data
%   This function generates synthetic CPU frequency data
%   for a Ryzen 5 5600x CPU (base: 3.7GHz, boost: 4.6GHz)
%
%   Parameters:
%       num_runs: Number of benchmark runs to generate
%
%   Each run is 1 minute long (at 10Hz sampling rate)
%   Failure types are evenly distributed
%   Data is saved to benchmark_data.csv

% CPU specifications
base_freq = 3.7; % Base frequency in GHz
boost_freq = 4.6; % Maximum boost frequency in GHz

% Sampling parameters
fs = 10; % Sampling frequency (Hz)
duration_seconds = 60;
num_samples = duration_seconds * fs;
all_data = [];

% Calculate number of runs for each failure type
runs_per_type = floor(num_runs / 4);
extra_runs = num_runs - (runs_per_type * 4);

failure_distribution = [];

% Add runs for each failure type
for failure_type = 0:3
    failure_distribution = [failure_distribution, ones(1, runs_per_type) * failure_type];
end

% Add any extra runs as no-failure
if extra_runs > 0
    failure_distribution = [failure_distribution, zeros(1, extra_runs)];
end

% Shuffle the failure distribution
rng('shuffle');  % Use different random seed each time
failure_distribution = failure_distribution(randperm(length(failure_distribution)));

%generate data for each run
for run = 1:num_runs
    t = (0:num_samples-1) / fs;   
    %get chosen failtype
    failure_type = failure_distribution(run);
    %generate frequency data
    freq_profile = generate_frequency_profile(num_samples, base_freq, boost_freq, fs);
    %inject fail
    if failure_type > 0
        freq_profile = inject_failure(freq_profile, num_samples, fs, failure_type);
    end
    %combine into single dataset
    run_data = [t(:), freq_profile(:), ones(num_samples, 1) * run, ones(num_samples, 1) * failure_type];
    %append
    all_data = [all_data; run_data];
end

% Create table and save to CSV
data_table = array2table(all_data, 'VariableNames', {'Time', 'Frequency_GHz', 'Run_ID', 'Failure_Type'});

% Convert numeric failure type to string description
failure_type_names = categorical(data_table.Failure_Type, 0:3, {'None', 'Thermal_Throttling', 'Frequency_Oscillation', 'Stuck_Frequency'});
data_table.Failure_Type = failure_type_names;

writetable(data_table, 'benchmark_data.csv');

% Count occurrences of each failure type
failure_counts = zeros(1, 4);
for i = 0:3
    failure_counts(i+1) = sum(all_data(:, 4) == i) / num_samples;
end

fprintf('Generated %d benchmark runs with distribution:\n', num_runs);
fprintf('  None: %d runs (%.1f%%)\n', failure_counts(1), failure_counts(1)/num_runs*100);
fprintf('  Thermal_Throttling: %d runs (%.1f%%)\n', failure_counts(2), failure_counts(2)/num_runs*100);
fprintf('  Frequency_Oscillation: %d runs (%.1f%%)\n', failure_counts(3), failure_counts(3)/num_runs*100);
fprintf('  Stuck_Frequency: %d runs (%.1f%%)\n', failure_counts(4), failure_counts(4)/num_runs*100);
fprintf('Data saved to benchmark_data.csv\n');

% Show example plot of the first run
run_ids = unique(all_data(:, 3));
first_run_data = all_data(all_data(:, 3) == run_ids(1), :);

figure;
plot(first_run_data(:, 1), first_run_data(:, 2));
title(sprintf('Run %d Frequency Profile (Failure: %s)', run_ids(1), char(failure_type_names(find(all_data(:, 3) == run_ids(1), 1)))));
xlabel('Time (s)');
ylabel('Frequency (GHz)');
ylim([base_freq-0.2, boost_freq+0.2]);
grid on;
end

%% Helper Function 1: Generate CPU Frequency Profile
function freq_profile = generate_frequency_profile(num_samples, base_freq, boost_freq, fs) 
    % Start around base frequency
    freq_profile = zeros(num_samples, 1);
    freq_profile(1) = base_freq + 0.1 * rand();
    
    % Define a standard benchmark pattern with some randomness
    % This simulates a consistent benchmark with 4 phases
    
    % Phase durations (in seconds)
    warmup_duration = 5;
    load1_duration = 15;
    mid_duration = 10;
    load2_duration = 25;
    cooldown_duration = 5;
    
    % Convert to samples
    warmup_samples = warmup_duration * fs;
    load1_samples = load1_duration * fs;
    mid_samples = mid_duration * fs;
    load2_samples = load2_duration * fs;
    cooldown_samples = cooldown_duration * fs;
    
    % Check if the sum matches expected duration (was added for debugging)
    total_samples = warmup_samples + load1_samples + mid_samples + load2_samples + cooldown_samples;
    if total_samples ~= num_samples
        cooldown_samples = cooldown_samples + (num_samples - total_samples);
    end
    
    % Generate target frequencies for each phase (with small random variation between runs)
    warmup_target = base_freq + 0.2 + 0.1 * rand();
    load1_target = boost_freq - 0.2 * rand();
    mid_target = base_freq + 0.5 + 0.2 * rand();
    load2_target = boost_freq - 0.1 * rand();
    cooldown_target = base_freq + 0.1 * rand();
    
    % Define phase boundaries
    phase1_end = warmup_samples;
    phase2_end = phase1_end + load1_samples;
    phase3_end = phase2_end + mid_samples;
    phase4_end = phase3_end + load2_samples;
    
    % Add some random variation to transitions (slight shifts in timing)
    phase_variation = round(0.5 * fs * rand(3, 1)); % Up to 0.5 second variation
    phase1_end = max(5, phase1_end + phase_variation(1));
    phase2_end = max(phase1_end + 5, phase2_end + phase_variation(2));
    phase3_end = max(phase2_end + 5, phase3_end + phase_variation(3));
    phase4_end = min(num_samples, phase4_end);
    
    % Calculate frequency changes with thermal throttling simulation
    freq_change_rate = 0.05; % Maximum GHz change per second
    max_step = freq_change_rate / fs; % Maximum change per sample
    
    % Fill in the frequency profile
    for i = 2:num_samples
        % Determine current target based on phase
        if i <= phase1_end
            target = warmup_target;
        elseif i <= phase2_end
            target = load1_target;
        elseif i <= phase3_end
            target = mid_target;
        elseif i <= phase4_end
            target = load2_target;
        else
            target = cooldown_target;
        end
        
        % Calculate direction and magnitude of change
        diff_to_target = target - freq_profile(i-1);
        step = sign(diff_to_target) * min(abs(diff_to_target), max_step);
        
        % Add small random variations (jitter)
        jitter = 0.02 * randn() / fs;
        
        freq_profile(i) = min(boost_freq, max(base_freq, freq_profile(i-1) + step + jitter));
    end
    
    % Add small high-frequency noise to simulate minor fluctuations
    freq_profile = freq_profile + 0.01 * randn(size(freq_profile));
    
    % Ensure within limits
    freq_profile = min(boost_freq, max(base_freq, freq_profile));
end

%% Helper Function 2: Inject Failures into CPU Data
function freq_profile = inject_failure(freq_profile, num_samples, fs, failure_type)
% Inject different types of failures into the data
    
    % Choose random failure start time (after 10% of the run)
    failure_start = randi([ceil(0.1 * num_samples), ceil(0.8 * num_samples)]);
    
    switch failure_type
        case 1
            % Type 1: Thermal throttling
            failure_duration = randi([2, 10]) * fs; % 2-10 seconds
            failure_end = min(failure_start + failure_duration, num_samples);
            
            % Drop frequency by 30-70%
            drop_factor = 0.3 + 0.4 * rand();
            freq_profile(failure_start:failure_end) = freq_profile(failure_start:failure_end) * drop_factor;
            
        case 2
            % Type 2: Frequency oscillation
            failure_duration = randi([5, 20]) * fs;
            failure_end = min(failure_start + failure_duration, num_samples);
            
            % Create oscillation pattern
            t = 0:(failure_end - failure_start);
            osc_freq = 0.5 + rand();
            osc_amp = 0.2 + 0.3 * rand();
            
            oscillation = osc_amp * sin(2 * pi * osc_freq * t / fs);
            freq_profile(failure_start:failure_end) = freq_profile(failure_start:failure_end) + oscillation';
            
            % Ensure within limits
            freq_profile = min(4.6, max(3.7, freq_profile));
            
        case 3
            % Type 3: Stuck frequency
            failure_duration = randi([10, 30]) * fs;
            failure_end = min(failure_start + failure_duration, num_samples);
            
            stuck_freq = freq_profile(failure_start) * (0.8 + 0.2 * rand());
            freq_profile(failure_start:failure_end) = stuck_freq;
    end
end