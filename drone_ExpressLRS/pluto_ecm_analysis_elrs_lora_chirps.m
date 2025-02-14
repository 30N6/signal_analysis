%unzip data before running

filenames = {"analysis-20250208-000401-500Hz.log", "500Hz";
             "analysis-20250208-000450-333Hz-full.log", "333Hz";
             "analysis-20250208-000540-250Hz.log", "250Hz";
             "analysis-20250208-000621-150Hz.log", "150Hz";
             "analysis-20250208-000703-100Hz.log", "100Hz";
             "analysis-20250208-000744-50Hz.log", "50Hz"};

Fs = 7.68e6;
L = 2048;

reload = 1;

if reload
    scan_reports = [];
    init_done = false;

    for i_filename = 1:length(filenames)
        filename = filenames{i_filename, 1};
        mode = filenames{i_filename, 2};

        fprintf("Reading %s...\n", filenames{i_filename});
        lines = readlines(filenames{i_filename});
        fprintf("Done reading %s.\n", filenames{i_filename});
        
        for ii = 1:length(lines)
            if strlength(lines(ii)) <= 1
                continue
            end

            if mod(ii, 1000) == 0
                fprintf("Processed %d lines\n", ii);
            end
        
            decoded_line = jsondecode(lines(ii));
            data = decoded_line.data;
            data.mode = mode;
            if ~isfield(data, 'iq_data')
                continue;
            end    
    
            if ~init_done
                scan_reports = data;
                init_done = true;
            else
                scan_reports(end + 1) = data;
            end
        end
    end
end

%%
for ii = 1:length(scan_reports)
    scan_reports(ii).iq_data = scan_reports(ii).iq_data(:, 1) + 1j * scan_reports(ii).iq_data(:, 2);
    scan_reports(ii).mean_power_dB = 20*log10(mean(abs(scan_reports(ii).iq_data)));
    scan_reports(ii).iq_data_padded = paddata(scan_reports(ii).iq_data, L);        

    scan_reports(ii).iq_phase = unwrap(atan2(imag(scan_reports(ii).iq_data), real(scan_reports(ii).iq_data)));
    scan_reports(ii).iq_freq = (1/(2*pi)) * diff(scan_reports(ii).iq_phase) / (1/Fs);        
    scan_reports(ii).timestamp_sec = scan_reports(ii).timestamp * (1/(4*61.44e6));

    if length(scan_reports(ii).iq_freq) < 50
        scan_reports(ii).r_squared_fsk_min = 0;
        scan_reports(ii).r_squared_fsk_max = 0;
    else
        l_freq = length(scan_reports(ii).iq_freq);
        rs = [0, 0];
        [rs(1), ~] = analyze_fsk(scan_reports(ii).iq_freq(1:floor(l_freq/2)), 2, Fs);
        [rs(2), ~] = analyze_fsk(scan_reports(ii).iq_freq(ceil(l_freq/2):end), 2, Fs);
        [~, fsk_len_hist, scan_reports(ii).fsk_freq_spread] = analyze_fsk(scan_reports(ii).iq_freq, 2, Fs);

        [~, scan_reports(ii).fsk_len_peak] = max(fsk_len_hist);
        scan_reports(ii).r_squared_fsk_min = min(rs);
        scan_reports(ii).r_squared_fsk_max = max(rs);
        
        [scan_reports(ii).fft_freq_mean, scan_reports(ii).fft_freq_std] = get_fft_stats(scan_reports(ii).iq_data_padded, Fs);

        [scan_reports(ii).lora_mean_slope, scan_reports(ii).lora_r_squared] = get_best_lora_poly_fit(scan_reports(ii).iq_freq, Fs);
        [scan_reports(ii).lora_peak_count_ratio, scan_reports(ii).lora_peak_spacing_ratio] = get_lora_fit_metrics(scan_reports(ii).iq_data_padded, scan_reports(ii).lora_mean_slope, Fs);
    end
end

scan_reports = scan_reports';

%% 
is_tx_listen = false(length(scan_reports), 1);
for ii = 1:length(scan_reports)
    is_tx_listen(ii) = scan_reports(ii).controller_state == "TX_LISTEN";
end
length_match    = ([scan_reports.iq_length] > 128).';
power_match     = ([scan_reports.mean_power_dB] > 10).';
timestamp_match = ([scan_reports.timestamp_sec] > 0).' & ([scan_reports.timestamp_sec] < 100).';
mod_match        = ([scan_reports.lora_r_squared] > 0.7).' & ([scan_reports.lora_peak_count_ratio] > 2).';
filtered_reports = scan_reports(length_match & is_tx_listen & power_match & timestamp_match & mod_match);

figure(1);
plot([filtered_reports.lora_mean_slope] * (1/(1e3/1e-6)), 'o');



% pulse_durations = unique([pdw_reports_filtered.pulse_duration]);
% pulse_duration_bins = min(pulse_durations):(max(pulse_durations)+1);
% pulse_duration_bin_x = pulse_duration_bins(1:end-1);
% pulse_duration_counts = zeros(length(pulse_frequencies_filtered), length(pulse_duration_bins)-1);
% pulse_duration_legend = strings(length(pulse_frequencies_filtered), 1);

modes = sort(unique([filtered_reports.mode]));

chirp_bins = -100:0.5:120;
chirp_bin_x = chirp_bins(1:end-1);
chirp_counts = zeros(length(modes), length(chirp_bins)-1);

for ii = 1:length(modes)
    mode = modes(ii);
    
    mode_match = false(length(filtered_reports), 1);
    for jj = 1:length(filtered_reports)
        mode_match(jj) = filtered_reports(jj).mode == mode;
    end
    data = filtered_reports(mode_match);
    current_chirps = [data.lora_mean_slope] * (1/(1e3/1e-6));
    
    chirp_counts(ii, :) = histcounts(current_chirps, chirp_bins);
    %pulse_duration_legend(ii) = sprintf('%0.2f MHz - median=%0.1f us', freq, median(current_pd) * pd_scale_factor_us);
end

figure(2);
bar(chirp_bin_x, chirp_counts.', 1, 'stacked');
grid on;
title('Chirp rate by ELRS mode');
xlabel('Chirp rate (kHz/us)');
ylabel('Count');
legend(modes);


%%
function [r_squared, det_length_hist, freq_spread] = analyze_fsk(data, M, Fs)
    freq_mean = mean(data);

    freq_cluster = ones([length(data), 1]);
    freq_cluster(data >= freq_mean) = 2;
    
    %freq_sorted = sort(data.iq_freq);
    %freq_span = freq_sorted(round(length(freq_sorted) * 0.9)) - freq_sorted(round(length(freq_sorted) * 0.1));
    %freq_transitions = (data.iq_freq < (freq_mean + 0.8*freq_span/2)) & (data.iq_freq > (freq_mean - 0.8*freq_span/2));
    %freq_cluster(freq_transitions) = 0;
    
    ss_res = [0, 0];
    cluster_y = [0, 0];
    
    det_length_hist = zeros([1, length(data)]);
    det_len = 1;
    for jj = 2:length(data)
        if freq_cluster(jj) == freq_cluster(jj - 1)
            det_len = det_len + 1;
        else
            det_length_hist(det_len) = det_length_hist(det_len) + 1;
            det_len = 1;
        end
    end

    for ii = 1:M
        y_k = data(freq_cluster == ii);
        cluster_y(ii) = mean(y_k); %median(y_k);
        ss_res(ii) = sum((y_k - cluster_y(ii)).^2);
    end
    y_tot = data(freq_cluster ~= 0);
    ss_tot = sum((y_tot - mean(y_tot)).^2);
    
    r_squared = 1 - sum(ss_res) / ss_tot;
    freq_spread = abs(diff(cluster_y));
end

function [freq_mean, freq_std] = get_fft_stats(data, Fs)
    L = length(data);
    f_shifted = (Fs/L)*(-L/2:L/2-1).';

    Y = fft(data);
    Y_shifted = fftshift(Y);
    
    Y_abs = abs(Y_shifted);
    Y_sum = sum(Y_abs);
    Y_scaled = Y_abs .* (1/Y_sum);
    Y_weighted_1 = Y_scaled .* f_shifted;

    freq_mean = sum(Y_weighted_1);
    Y_var = sum(Y_scaled .* (f_shifted - freq_mean).^2);

    freq_std = sqrt(Y_var);
end

function [mean_slope, r_squared] = try_lora_poly_fit(iq_freq, num_chunks, Fs)
    chunk_length = floor(length(iq_freq)/num_chunks);

    freq_chunks = zeros([num_chunks, chunk_length]);
    p_freq      = zeros([num_chunks, chunk_length]);
    r_squared   = zeros([num_chunks, 1]);
    freq_poly   = zeros([num_chunks, 2]);

    for ii = 1:num_chunks
        freq_chunks(ii, :) = iq_freq((ii-1)*chunk_length+1 : ii*chunk_length);

        x_freq = (0:chunk_length-1) * (1/Fs);
        [freq_poly(ii, :), S_freq] = polyfit(x_freq, freq_chunks(ii, :), 1);
        ss_res = S_freq.normr ^ 2;
        ss_tot = sum((freq_chunks(ii, :) - mean(freq_chunks(ii, :))).^2);
        r_squared(ii) = 1 - ss_res/ss_tot;
        p_freq(ii, :) = freq_poly(ii, 1) * x_freq + freq_poly(ii, 2);        
    end
    
    mean_slope = freq_poly(:, 1);
end

function [mean_slope, r_squared] = get_best_lora_poly_fit(iq_freq, Fs)
    num_chunks = [4, 16];
    mean_slope = [];
    r_squared = [];

    for ii = 1:length(num_chunks)
        [ms, rs] = try_lora_poly_fit(iq_freq, num_chunks(ii), Fs);
        mean_slope = [mean_slope; ms];
        r_squared = [r_squared; rs];
    end
    [~, max_i] = max(r_squared);
    r_squared = r_squared(max_i);
    mean_slope = mean_slope(max_i);
end

function [peak_count_ratio, peak_spacing_ratio] = get_lora_fit_metrics(iq_data, original_slope, Fs)
    peak_threshold = 0.25;

    L = length(iq_data);

    Y_original = fftshift(abs(fft(iq_data)));

    t = (0:L-1).' * (1/Fs);
    f_chirp = -0.5 * original_slope * t;
    y_chirp = exp(1j*2*pi*f_chirp.*t);
    
    iq_dechirped = iq_data .* y_chirp;
    iq_dechirped = paddata(iq_dechirped, L);

    Y_dechirped = fftshift(abs(fft(iq_dechirped)));

    f_a = find(Y_original > peak_threshold * max(Y_original));
    y_a = Y_original(Y_original > peak_threshold * max(Y_original));
    f_b = find(Y_dechirped > peak_threshold * max(Y_dechirped));
    y_b = Y_dechirped(Y_dechirped > peak_threshold * max(Y_dechirped));
    
    if (length(y_a) > 2)
        [peak_a, loc_a] = findpeaks(y_a, f_a);
    else
        [peak_a, loc_a] = max(y_a);
    end

    if (length(y_b) > 2)
        [peak_b, loc_b] = findpeaks(y_b, f_b);
    else
        [peak_b, loc_b] = max(y_b);
    end
    
    peak_count_ratio = length(peak_a)/length(peak_b);
    if length(loc_b) > 1
        peak_spacing_ratio = mean(diff(loc_a)) / mean(diff(loc_b));
    else
        peak_spacing_ratio = 0;
    end
end