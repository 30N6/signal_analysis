%filename = 'analysis-20250118-164133-10mW-near-50ms.log'; pri_raster_period_us = 41666.954;
%filename = 'analysis-20250118-171315-250mW-near-50ms-50Hz.log'; pri_raster_period_us = 41666.954;

filename = 'analysis-20250118-165648-250mW-near-50ms-F1000.log'; pri_raster_period_us = 1000;
%filename = 'analysis-20250118-173017-250mW-near-50ms-50Hz.log'; pri_raster_period_us = 40e3;
%filename = 'analysis-20250118-121543-10mW-far.log'; pri_raster_period_us = 1000;

min_pulses_per_freq         = 0;
min_pulse_power_dB          = 10;
min_duration_us             = 1; %50;
require_modulation          = true;

toa_scale_factor_s          = 1/(4*61.44e6);
pd_scale_factor_us          = 1/1.92;
channel_spacing             = 0.96;
max_toa_diff_levels         = 200;
max_analysis_pri_single_us  = 20e3;
max_analysis_pri_full_us    = 100e3;

pri_num_clusters            = 16;

pri_raster_offset_us        = 1000;
pri_raster_num_clusters     = 12;
pri_cluster_min_pulses      = 0.05;

reload = 1;
if reload
    lines = readlines(filename);
    pdw_reports = [];
    init_done = false;
    
    for ii = 1:length(lines)
        if strlength(lines(ii)) <= 1
            continue
        end
    
        decoded_line = jsondecode(lines(ii));
        report = decoded_line.data;
        if ~isfield(report, 'pulse_seq_num')
            continue;
        end    
        
        if ~report.buffered_frame_valid
            report.buffered_frame_data = zeros([50, 2]);
        end

        if ~init_done
            pdw_reports = report;
            init_done = true;
        else
            pdw_reports(end + 1) = report;
        end
    end
    clear('lines');

    pdw_reports = pdw_reports';

    if (require_modulation)
        filtered_reports = [];

        for ii = 1:length(pdw_reports)
            if isstruct(pdw_reports(ii).modulation_data) && strcmp(pdw_reports(ii).modulation_data.modulation_type, 'FM')
                if isempty(filtered_reports)
                    filtered_reports = pdw_reports(ii);
                else
                    filtered_reports(end + 1) = pdw_reports(ii);
                end
            end
        end

        pdw_reports = filtered_reports.';
        clear('filtered_reports');
    end
    
    toa_0 = pdw_reports(1).pulse_start_time * toa_scale_factor_s;

    for ii = 1:length(pdw_reports)
        pdw_reports(ii).pulse_toa_s             = pdw_reports(ii).pulse_start_time * toa_scale_factor_s - toa_0;
        pdw_reports(ii).pulse_duration_us       = pdw_reports(ii).pulse_duration * pd_scale_factor_us;

        pdw_reports(ii).recorded_iq_data        = 1j * pdw_reports(ii).buffered_frame_data(:, 1) + pdw_reports(ii).buffered_frame_data(:, 2);
        pdw_reports(ii).recorded_power          = abs(pdw_reports(ii).recorded_iq_data) .^ 2;
        pdw_reports(ii).recorded_noise_power    = mean(pdw_reports(ii).recorded_power(1:8));
    
        iq_length = min(50 - 8, pdw_reports(ii).pulse_duration);
        pdw_reports(ii).recorded_pulse_power    = mean(pdw_reports(ii).recorded_power(9:(9 + iq_length - 1)));
        pdw_reports(ii).recorded_pulse_snr      = pdw_reports(ii).recorded_pulse_power / max(pdw_reports(ii).recorded_noise_power, 1);
        pdw_reports(ii).pulse_power             = pdw_reports(ii).pulse_power_accum / pdw_reports(ii).pulse_duration;
        pdw_reports(ii).implied_pulse_snr       = pdw_reports(ii).pulse_power / pdw_reports(ii).pulse_threshold;
    end

    original_reports = pdw_reports;
end

valid_reports = (10*log10([original_reports.pulse_power]) > min_pulse_power_dB) & ([original_reports.pulse_duration_us] > min_duration_us);
pdw_reports = original_reports(valid_reports);

[~, sorted_i] = sort([pdw_reports.pulse_start_time]);
pdw_reports = pdw_reports(sorted_i);

pulse_frequencies = sort(unique([pdw_reports.channel_frequency]));
pdw_reports_by_freq = dictionary;
for ii = 1:length(pulse_frequencies)
    freq = pulse_frequencies(ii);
    valid_reports = ([pdw_reports.channel_frequency] == freq);
    if sum(valid_reports) < min_pulses_per_freq
        continue
    end
    pdw_reports_by_freq{freq} = pdw_reports(valid_reports);
end

pulse_frequencies_filtered = pdw_reports_by_freq.keys;
pdw_reports_filtered = pdw_reports(ismember([pdw_reports.channel_frequency], pulse_frequencies_filtered));

freq_legend = strings(length(pulse_frequencies_filtered), 1);
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    freq_legend(ii) = sprintf('%0.2f MHz', freq);
end

freq_hist_bins = min(pulse_frequencies) : channel_spacing : (max(pulse_frequencies) + channel_spacing);

figure(1);
histogram([pdw_reports.channel_frequency], freq_hist_bins);
grid on;
title('Pulse frequency');
xlabel('Frequency (MHz)');
ylabel('Count');


pulse_durations = unique([pdw_reports_filtered.pulse_duration]);
pulse_duration_bins = min(pulse_durations):(max(pulse_durations)+1);
pulse_duration_bin_x = pulse_duration_bins(1:end-1);
pulse_duration_counts = zeros(length(pulse_frequencies_filtered), length(pulse_duration_bins)-1);
pulse_duration_legend = strings(length(pulse_frequencies_filtered), 1);
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    current_pd = [pdw_reports_by_freq{freq}.pulse_duration];
    pulse_duration_counts(ii, :) = histcounts(current_pd, pulse_duration_bins);
    pulse_duration_legend(ii) = sprintf('%0.2f MHz - median=%0.1f us', freq, median(current_pd) * pd_scale_factor_us);
end

figure(2);
bar(pulse_duration_bin_x * pd_scale_factor_us, pulse_duration_counts.', 1, 'stacked');
grid on;
title('Pulse duration by frequency');
xlabel('Duration (us)');
ylabel('Count');
legend(pulse_duration_legend);


figure(3);
hold off;
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    reports = pdw_reports_by_freq{freq};
    power_dB = 10*log10([reports.pulse_power]);
    plot([reports.pulse_toa_s], power_dB, '.-');
    hold on;
    %pulse_toa       = [pdw_reports_by_freq{freq}.pulse_toa_us];
    %pulse_duration  = [pdw_reports_by_freq{freq}.pulse_duration_us];
end
grid on;
title('Pulse power');
xlabel('Time (s)');
ylabel('Power (dB)');
legend(freq_legend);


pulse_toa_us = [pdw_reports.pulse_toa_s] / 1e-6;
pulse_dwell_seq = [pdw_reports.dwell_seq_num];
dwell_seqs = unique(pulse_dwell_seq);

pulse_toa_diff_full = [];
pulse_toa_diff_single = [];
for jj = 1:length(dwell_seqs)
    dwell_seq = dwell_seqs(jj);
    matching_toa = pulse_toa_us(pulse_dwell_seq == dwell_seq);

    toa_diff = get_toa_diff_single_dwell(matching_toa, max_toa_diff_levels);
    pulse_toa_diff_full = [pulse_toa_diff_full; toa_diff];
    toa_diff = get_toa_diff_single_dwell(matching_toa, 1);
    pulse_toa_diff_single = [pulse_toa_diff_single; toa_diff];        
end

filtered_pulse_toa_diff_full = pulse_toa_diff_full(pulse_toa_diff_full < max_analysis_pri_full_us);
filtered_pulse_toa_diff_single = pulse_toa_diff_single(pulse_toa_diff_single < max_analysis_pri_single_us);

figure(4);
histogram(filtered_pulse_toa_diff_full, 1000);
grid on;
xlabel('PRI (us)');
ylabel('Count');
title(sprintf('Full PRI histogram: %0.2f MHz', freq));

figure(5);
histogram(filtered_pulse_toa_diff_single, 1000);    
grid on;
xlabel('PRI (us)');
ylabel('Count');
title(sprintf('First-level PRI histogram: %0.2f MHz', freq));    

% figure(6);
% ax(ii, 3) = subplot(length(pulse_frequencies_filtered), 1, ii);
% toa_all = [reports.pulse_toa_s] / 1e-6;
% toa_all_diff = diff(toa_all);
% plot(toa_all(1:end-1), toa_all_diff, 'o');

%% PRI raster

figure(10);
pulse_toa_us = [pdw_reports.pulse_toa_s] / 1e-6 + pri_raster_offset_us;

toa_x = mod(pulse_toa_us, pri_raster_period_us);
toa_y = floor(pulse_toa_us / pri_raster_period_us) * pri_raster_period_us * 1e-6;
plot(toa_x, toa_y, 'o');
grid on;
xlabel('Fast time (us)');
ylabel('Slow time (s)');
title(sprintf('PRI raster: period=%0.2f us', pri_raster_period_us));
legend(freq_legend);

figure(11);  
subplot(2,1,1);
pulse_toa_us = [pdw_reports.pulse_toa_s] / 1e-6 + pri_raster_offset_us;
toa_x = mod(pulse_toa_us, pri_raster_period_us).';
toa_y = floor(pulse_toa_us / pri_raster_period_us).' * pri_raster_period_us * 1e-6;
cluster_idx = kmeans(toa_x, pri_raster_num_clusters);

cluster_median = zeros(pri_raster_num_clusters, 1);
cluster_count = zeros(pri_raster_num_clusters, 1);
for ii = 1:pri_raster_num_clusters
    current_cluster_valid = cluster_idx == ii;
    cluster_y = toa_y(current_cluster_valid);
    cluster_x = toa_x(current_cluster_valid);
    cluster_median(ii) = median(cluster_x);
    cluster_count(ii) = length(cluster_x);

    plot(cluster_x, cluster_y, 'o', [cluster_median(ii), cluster_median(ii)], [cluster_y(1), cluster_y(end) + 2], '-');
    hold on;
end

ax1 = subplot(2,1,1);
grid on;
xlabel('Fast time (us)');
ylabel('Slow time (s)');
title(sprintf('PRI raster clustering: period=%0.2f us', pri_raster_period_us));

ax2 = subplot(2,1,2);
stem(cluster_median, cluster_count);
grid on;
xlabel('Fast time (us)');
ylabel('Count');
title('PRI raster cluster counts');
linkaxes([ax1, ax2], 'x');

fprintf('Stagger pattern from raster, freq=%.2f:\n', pri_analysis_freq);
sorted_median = sort(cluster_median);
pri_stagger = zeros(length(sorted_median), 1);
for ii = 1:pri_raster_num_clusters
    i_0 = ii;
    i_1 = mod(ii, pri_raster_num_clusters) + 1;

    current_pri = mod(sorted_median(i_1) - sorted_median(i_0), pri_raster_period_us);
    pri_stagger(ii) = current_pri;

    fprintf('  %2d: PRI: %0.2f us   PRF: %0.2f Hz\n', ii, current_pri, 1/(current_pri * 1e-6));
end

%%
return
disp('Pausing for figure adjustments');
pause();

%% 
figs_to_save = [1, 2, 3, 4, 5, 9, 10, 11];
split_fn = split(filename, '.');
filename_base = split_fn{1};

for ii = 1:length(figs_to_save)
    fig_index = figs_to_save(ii);
    f = figure(fig_index);
    f.Position = [1400 100 800 500];
    
    fig_filename = sprintf('%s_fig_%d.png', filename_base, fig_index);
    saveas(f, fig_filename);
end

%%
function td = get_toa_diff_single_dwell(toa, N_diffs)
    td = [];

    for i_diff = 1:N_diffs
        current_td = zeros((length(toa) - i_diff), 1);
        for ii = 1:(length(toa) - i_diff)
            current_td(ii) = toa(ii + i_diff) - toa(ii);
        end
        td = [td; current_td];
    end
end
