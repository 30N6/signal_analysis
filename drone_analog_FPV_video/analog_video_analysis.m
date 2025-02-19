%filename = 'analysis-20250118-135026-5917-20mW-250ms.log'; pri_raster_period_us = 33365.4; min_pulses_per_freq = 5000;
filename = 'analysis-20250118-143522-5917-20mW-10ms.log'; pri_raster_period_us = 33365.8; min_pulses_per_freq = 2000;

toa_scale_factor_s          = 1/(4*61.44e6);
pd_scale_factor_us          = 1/1.92;
channel_spacing             = 0.96;
max_toa_diff_levels         = 12;
max_analysis_pri_full_us    = 50e3;
max_analysis_pri_single_us  = 500;

pri_num_clusters            = 16;

pri_raster_offset_us        = 1000;
pri_raster_num_clusters     = 12;
pri_cluster_min_pulses      = 0.05;

%pri_analysis_freq           = 5917.44;
pri_analysis_freq           = 5918.4;
%pri_analysis_freq           = 1056.96;

filter_freq_range = [5900, 6000];

reload = 0;
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

    pdw_reports = pdw_reports';
    
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
    
    pulse_frequencies = sort(unique([pdw_reports.channel_frequency]));
    
    pdw_reports_by_freq = dictionary;
    for ii = 1:length(pulse_frequencies)
        freq = pulse_frequencies(ii);
        valid_reports = [pdw_reports.channel_frequency] == freq;
        if (freq < filter_freq_range(1)) || (freq > filter_freq_range(2))
            continue
        end
        if sum(valid_reports) < min_pulses_per_freq
            continue
        end
        pdw_reports_by_freq{freq} = pdw_reports(valid_reports);
    end
    
    pulse_frequencies_filtered = pdw_reports_by_freq.keys;
    pdw_reports_filtered = pdw_reports(ismember([pdw_reports.channel_frequency], pulse_frequencies_filtered));
end

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
power_legend = strings(length(pulse_frequencies_filtered), 1);
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    reports = pdw_reports_by_freq{freq};
    power_dB = 10*log10([reports.pulse_power]);
    mean_power_dB = 10*log10(mean([reports.pulse_power]));
    power_legend(ii) = sprintf('%0.2f MHz - mean=%0.1f dB', freq, mean_power_dB);

    plot([reports.pulse_toa_s], power_dB, '.-');
    hold on;
    %pulse_toa       = [pdw_reports_by_freq{freq}.pulse_toa_us];
    %pulse_duration  = [pdw_reports_by_freq{freq}.pulse_duration_us];
end
grid on;
title('Pulse power');
xlabel('Time (s)');
ylabel('Power (dB)');
legend(power_legend);


ax = zeros(length(pulse_frequencies_filtered), 3);
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    reports = pdw_reports_by_freq{freq};

    pulse_toa_us = [reports.pulse_toa_s] / 1e-6;
    pulse_dwell_seq = [reports.dwell_seq_num];
    dwell_seqs = unique(pulse_dwell_seq);

    pulse_toa_diff_single = [];
    for jj = 1:length(dwell_seqs)
        dwell_seq = dwell_seqs(jj);
        matching_toa = pulse_toa_us(pulse_dwell_seq == dwell_seq);

        toa_diff = get_toa_diff_single_dwell(matching_toa, 1);
        pulse_toa_diff_single = [pulse_toa_diff_single; toa_diff];        
    end
    
    filtered_pulse_toa_diff_single = pulse_toa_diff_single(pulse_toa_diff_single < max_analysis_pri_single_us);
    
    figure(5);
    ax(ii, 2) = subplot(length(pulse_frequencies_filtered), 1, ii);
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
end
linkaxes(ax(:, 2), 'x');
% linkaxes(ax(:, 3), 'xy');

%% PRI cluster analysis

ax = zeros(length(pulse_frequencies_filtered), 3);
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    reports = pdw_reports_by_freq{freq};

    pulse_toa_us = [reports.pulse_toa_s] / 1e-6;
    pulse_toa_diff = diff(pulse_toa_us);
    pulse_toa_diff_valid = pulse_toa_diff < max_analysis_pri_single_us;
    
    pulse_toa_filt = pulse_toa_us([pulse_toa_diff_valid, false]).';
    pulse_toa_diff_filt = pulse_toa_diff(pulse_toa_diff_valid).';

    cluster_idx = kmeans(pulse_toa_diff_filt, pri_num_clusters);
    
    cluster_median = zeros(pri_num_clusters, 1);
    cluster_count = zeros(pri_num_clusters, 1);
    for jj = 1:pri_num_clusters
        current_cluster_valid = cluster_idx == jj;
        current_toa = pulse_toa_filt(current_cluster_valid);
        current_diff = pulse_toa_diff_filt(current_cluster_valid);
        cluster_median(jj) = median(current_diff);
        cluster_count(jj) = length(current_diff);

        if freq == pri_analysis_freq
            figure(9);  
            subplot(2,1,1);
            plot(current_diff, current_toa * 1e-6, 'o', [cluster_median(jj), cluster_median(jj)], [current_toa(1) * 1e-6, current_toa(end) * 1e-6 + 2], '-');
            hold on;
        end
    end
    
    cluster_threshold = pri_cluster_min_pulses * length(pulse_toa_filt);
    
    if freq == pri_analysis_freq
        figure(9);  
        ax1 = subplot(2,1,1);
        grid on;
        title('PRI clusters vs time');
        xlabel('PRI (us)');
        ylabel('Time (s)');

        ax2 = subplot(2,1,2);
        hold off;
        stem(cluster_median, cluster_count);
        hold on;
        plot([min(cluster_median), max(cluster_median)], [cluster_threshold, cluster_threshold], '--');
        grid on;
        title('PRI cluster pulse count');
        xlabel('PRI (us)');
        ylabel('Count');
        legend('PRI count', sprintf('Threshold: %0.1f%%', pri_cluster_min_pulses * 100));
        linkaxes([ax1, ax2], 'x');
    end
end

%% PRI raster

figure(10);
hold off;
for ii = 1:length(pulse_frequencies_filtered)
    freq = pulse_frequencies_filtered(ii);
    reports = pdw_reports_by_freq{freq};
    pulse_toa_us = [reports.pulse_toa_s] / 1e-6 + pri_raster_offset_us;
    
    toa_x = mod(pulse_toa_us, pri_raster_period_us);
    toa_y = floor(pulse_toa_us / pri_raster_period_us) * pri_raster_period_us * 1e-6;
    plot(toa_x, toa_y, 'o');
    hold on;
end
grid on;
xlabel('Fast time (us)');
ylabel('Slow time (s)');
title(sprintf('PRI raster: period=%0.2f us', pri_raster_period_us));
legend(freq_legend);


figure(11);
reports = pdw_reports_by_freq{pri_analysis_freq};
pulse_toa_us = [reports.pulse_toa_s] / 1e-6 + pri_raster_offset_us;

toa_x = mod(pulse_toa_us, pri_raster_period_us);
toa_y = floor(pulse_toa_us / pri_raster_period_us) * pri_raster_period_us * 1e-6;
plot(toa_x, toa_y, 'o');
grid on;
xlabel('Fast time (us)');
ylabel('Slow time (s)');
title(sprintf('PRI raster: freq=%0.2f period=%0.2f us', pri_analysis_freq, pri_raster_period_us));


%%
figs_to_save = [1, 2, 3, 5, 10, 11];
for ii = 1:length(figs_to_save)
    fig_index = figs_to_save(ii);
    f = figure(fig_index);
    f.Position = [1400 100 800 500];
end

disp('Pausing for figure adjustments');
pause();

%% 

split_fn = split(filename, '.');
filename_base = split_fn{1};

for ii = 1:length(figs_to_save)
    fig_index = figs_to_save(ii);
    f = figure(fig_index);
    
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
