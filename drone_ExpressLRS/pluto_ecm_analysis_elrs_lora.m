%filename = "analysis-20250208-000007-F1000.log";
%filename = "analysis-20250208-000110-F500.log";
%filename = "analysis-20250208-000211-D500.log";
%filename = "analysis-20250208-000307-D250.log";
filename = "analysis-20250208-000401-500Hz.log";
%filename = "analysis-20250208-000450-333Hz-full.log";
%filename = "analysis-20250208-000540-250Hz.log";
%filename = "analysis-20250208-000621-150Hz.log";
%filename = "analysis-20250208-000703-100Hz.log";
%filename = "analysis-20250208-000744-50Hz.log";

%filename = "analysis-20250208-005119.log"; %background
%filename = "analysis-20250208-105833-F1000-10mW.log";

Fs = 7.68e6;
dt = 1/Fs;
L = 2048;

reload = 1;

if reload
    lines = readlines(filename);
    scan_reports = [];
    init_done = false;
    
    for ii = 1:length(lines)
        if strlength(lines(ii)) <= 1
            continue
        end
    
        decoded_line = jsondecode(lines(ii));
        data = decoded_line.data;
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

    for ii = 1:length(scan_reports)
        scan_reports(ii).iq_data = scan_reports(ii).iq_data(:, 1) + 1j * scan_reports(ii).iq_data(:, 2);
        %scan_reports(ii).iq_data = scan_reports(ii).iq_data - mean(scan_reports(ii).iq_data) * 1.0; %0.999;
        %if (mod(scan_reports(ii).channel_index, 2) == 1)
        %    scan_reports(ii).iq_data = baseband(scan_reports(ii).iq_data);
        %end
        scan_reports(ii).mean_power_dB = 20*log10(mean(abs(scan_reports(ii).iq_data)));
        scan_reports(ii).iq_data_padded = paddata(scan_reports(ii).iq_data, L);        

        scan_reports(ii).iq_phase = unwrap(atan2(imag(scan_reports(ii).iq_data), real(scan_reports(ii).iq_data)));
        scan_reports(ii).iq_freq = (1/(2*pi)) * diff(scan_reports(ii).iq_phase) / dt;        
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
end


freqs = unique([scan_reports.dwell_freq]);
figure(10);
for ii = 1:length(freqs)
    subplot(length(freqs), 1, ii)
    r = scan_reports([scan_reports.dwell_freq] == freqs(ii));
    plot([r.timestamp] * (1/(4*61.44e6)), [r.mean_power_dB], 'o');
end

freqs = unique([scan_reports.dwell_freq]);
figure(11);
for ii = 1:length(freqs)
   subplot(length(freqs), 1, ii)
   r = scan_reports([scan_reports.dwell_freq] == freqs(ii));
   plot([r.timestamp] * (1/(4*61.44e6)), [r.r_squared_fsk_min], 'o', [r.timestamp] * (1/(4*61.44e6)), [r.r_squared_fsk_max], 'o');
end



filter_freq = 2475;

is_tx_listen = false(length(scan_reports), 1);
for ii = 1:length(scan_reports)
    is_tx_listen(ii) = scan_reports(ii).controller_state == "TX_LISTEN";
end
freq_match      = ([scan_reports.dwell_freq] == filter_freq).';
length_match    = ([scan_reports.iq_length] > 128).';
power_match     = ([scan_reports.mean_power_dB] > 10).';
timestamp_match = ([scan_reports.timestamp_sec] > 0).' & ([scan_reports.timestamp_sec] < 100).';
mod_match        = ([scan_reports.lora_r_squared] > 0.7).' & ([scan_reports.lora_peak_count_ratio] > 2).';
filtered_reports = scan_reports(freq_match & length_match & is_tx_listen & power_match & timestamp_match & mod_match);

figure(20);
subplot(3,1,1);
%plot([filtered_reports.mean_power_dB]);
%plot([filtered_reports.timestamp] * (1/(4*61.44e6)), [filtered_reports.mean_power_dB], 'o');
plot([filtered_reports.mean_power_dB], 'o');
subplot(3,1,2);
plot([filtered_reports.iq_length], [filtered_reports.mean_power_dB], 'o');

figure(21);
plot([filtered_reports.lora_mean_slope] * (1/(1e3/1e-6)), 'o');

%%
offset = 0;
num_rows = 6;
num_cols = 6;

figure(2);
ax1 = zeros(num_rows, num_cols);
figure(3);
ax2 = zeros(num_rows, num_cols);

for row = 1:num_rows
    for col = 1:num_cols
        plot_index = (row-1) * num_cols + col;
        d = filtered_reports(plot_index + offset);        

        channel_freq = d.dwell_freq + (Fs/2)/1e6 * (d.channel_index - 8);
        %s = sprintf("[%d] %d %d: %.1f %0.3f %0.3f", plot_index + offset, d.dwell_freq, d.channel_index, channel_freq, d.r_squared_fsk_max, d.r_squared_fsk_min);
        s = sprintf("[%d]: %.1f %0.3f %0.3f %0.3f", plot_index + offset, channel_freq, d.lora_r_squared, d.lora_peak_count_ratio, d.lora_peak_spacing_ratio);
        
        t = (0:d.iq_length-1).' / Fs;
        y = d.iq_data(1:d.iq_length);
        
        figure(2);
        ax1(row, col) = subplot(num_rows,num_cols,plot_index);
        plot(t, real(y), t, imag(y));
        title(s);
    end
end

for row = 1:num_rows
    for col = 1:num_cols
        plot_index = (row-1) * num_cols + col;
        d = filtered_reports(plot_index + offset);        

        channel_freq = d.dwell_freq + (Fs/2)/1e6 * (d.channel_index - 8);
        %s = sprintf("[%d] %d %d: %.1f %0.3f %0.3f", plot_index + offset, d.dwell_freq, d.channel_index, channel_freq, d.r_squared_fsk_max, d.r_squared_fsk_min);
        s = sprintf("[%d]: %.1f %0.3f %0.3f %0.3f", plot_index + offset, channel_freq, d.lora_r_squared, d.lora_peak_count_ratio, d.lora_peak_spacing_ratio);
        
        t = (0:d.iq_length-1).' / Fs;
        y = d.iq_data(1:d.iq_length);
    
        figure(3);
        ax2(row, col) = subplot(num_rows,num_cols,plot_index);
        plot(t(1:end-1), d.iq_freq);
        title(s);
    end
end

%%
figure(30);
hold off;
for ii = 1:36
    d = filtered_reports(ii + offset);  
    t = (0:d.iq_length-1).' / Fs;
    y = d.iq_data(1:d.iq_length);    
    plot(t(1:end-1), d.iq_freq,'.');
    hold on;
end
%linkaxes(ax, 'y');


%%
d = filtered_reports(2);

rows = 6;

t = (0:d.iq_length-1).' / Fs;
y = d.iq_data(1:d.iq_length);

t_padded = (0:L-1).' / Fs;
y_padded = d.iq_data_padded;

figure(4);
ax_1 = subplot(rows,2,1);
plot(t, real(y), t, imag(y));

ax_2 = subplot(rows,2,2);
instfreq(y, Fs);
%linkaxes([ax_1, ax_2], 'x');

subplot(rows,2,3);
[c,lags] = xcorr(y);
plot(lags * (1/Fs) * 1e6, (abs(c)));

subplot(rows,2,4); 
X = fft(y_padded);
xc = abs(ifft(X .* conj(X)));

tx = t_padded(1:L/2);
xc = xc(1:L/2);

m_xc = mean(xc);
s_xc = std(xc) * 1.0;

i_th = xc > (m_xc + s_xc);
xc_th = xc(i_th);
t_th = tx(i_th);

plot(tx, xc, [0, tx(end)], [m_xc, m_xc], [0, tx(end)], [m_xc + s_xc, m_xc + s_xc], t_th, xc_th, 'o');

subplot(rows,2,5);
sf = compute_sfft(y_padded, 16);
imagesc(sf.');

subplot(rows,2,6);
Y = fft(y);
%plot(Fs/L*(-L/2:L/2-1), 20*log10(abs(fftshift(Y))),"LineWidth",1)
%plot(Fs/L*(-L/2:L/2-1), (abs(fftshift(Y))),"LineWidth",1)
plot(abs(fftshift(Y)));

subplot(rows,2,7);
plot(t, d.iq_phase);

subplot(rows,2,8);
plot(t(1:end-1), d.iq_freq);



freq_diff = diff(d.iq_freq);
freq_diff_clipped = freq_diff;
freq_diff_clipped(freq_diff_clipped > 1e5) = 0;
freq_diff_clipped(freq_diff_clipped < -1e5) = 0;

subplot(rows,2,9);
plot(t(1:end-2), freq_diff)

subplot(rows,2,10);
plot(t(1:end-2), freq_diff_clipped);

subplot(rows,2,11);
plot(t(1:end-2), cumsum(freq_diff_clipped));

subplot(rows, 2, 12);
plot(filter(ones(16,1), 1, d.iq_freq));

fprintf("corr ratio: %f\n", s_xc/m_xc);

%%

d = filtered_reports(8);
%d = filtered_reports(10);

t_pad = (0:L-1).' / Fs;
y_pad = d.iq_data_padded;

t_freq = (0:length(d.iq_freq)-1).' * (1/Fs);
y_freq = d.iq_freq;

%[r_squared, det_length_hist, freq_spread_fsk] = analyze_fsk(d.iq_freq, 2, Fs);
[fft_mean, fft_std] = get_fft_stats(d.iq_data_padded, Fs);
[mean_slope, iq_dechirped, f_dechirped] = analyze_lora(d.iq_data, d.iq_freq, Fs);

[best_mean_slope, best_r_squared] = get_best_lora_poly_fit(d.iq_freq, Fs);
[peak_count_ratio, peak_spacing_ratio] = get_lora_fit_metrics(d.iq_data_padded, best_mean_slope, Fs);

Y = fft(y_pad);
Y_shifted = abs(fftshift(Y));
f_shifted = (Fs/L)*(-L/2:L/2-1).';

figure(5);
subplot(4,2,1);
plot(t_pad, real(y_pad), t_pad, imag(y_pad));

subplot(4,2,2);
plot(t_freq, y_freq);

ax1 = subplot(4, 2, 3);
y_max = max(Y_shifted);
plot(f_shifted, Y_shifted, [fft_mean, fft_mean], [0, y_max], ...
    [fft_mean - fft_std, fft_mean - fft_std], [0, y_max], ...
    [fft_mean + fft_std, fft_mean + fft_std], [0, y_max]);
%TODO: mean, std
%f_shifted, cumsum(Y_scaled)

circ_xcorr = abs(ifft(Y .* conj(Y)));
subplot(4, 2, 4);
plot(circ_xcorr(1:L/2));

ax2 = subplot(4, 2, 5);
iq_dechirped_padded = paddata(iq_dechirped, L);
Y_dechirped = abs(fft(iq_dechirped_padded));
plot(f_shifted, fftshift(Y_dechirped));

subplot(4, 2, 6);
plot(t_freq, iq_dechirped);

ax3 = subplot(4,2,7);
plot(f_shifted(1:end-1), diff(Y_shifted), f_shifted(1:end-1), diff(fftshift(Y_dechirped)));
 
linkaxes([ax1, ax2, ax3], "x");

subplot(4,2,8);
thresh = 0.25;
f_a = find(Y_shifted > thresh * max(Y_shifted));
y_a = Y_shifted(Y_shifted > thresh * max(Y_shifted));
f_b = find(Y_dechirped > thresh * max(Y_dechirped));
y_b = Y_dechirped(Y_dechirped > thresh * max(Y_dechirped));
[peak_a, loc_a] = findpeaks(y_a, f_a);
[peak_b, loc_b] = findpeaks(y_b, f_b);

fprintf("peak_count: %d/%d=%0.2f\n", length(peak_a), length(peak_b), length(peak_a)/length(peak_b));
fprintf("loc_spacing: %0.1f/%0.1f=%0.2f\n", mean(diff(loc_a)), mean(diff(loc_b)), mean(diff(loc_a)) / mean(diff(loc_b)));

plot(f_a, y_a, 'o', loc_a, peak_a, 'x');
hold on;
plot(f_b, y_b, 'o', loc_b, peak_b, 'x');
hold off;

%plot(1:length(Y_shifted), cumsum(sort(Y_shifted)), 1:length(Y_dechirped), cumsum(sort(Y_dechirped)));
%plot()



%%
function r = compute_sfft(data, N)
    r = zeros(length(data)/N, N);
    for ii = 1:(length(data)/N)
        r(ii, :) = fftshift(abs(fft(data((ii-1)*N+1 : ii*N))));
    end
end

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

function [mean_slope, iq_dechirped, f_dechirped] = analyze_lora(iq_data, iq_freq, Fs)
    num_chunks = 4;

    r_squared_threshold = 0.9;

    chunk_length = floor(length(iq_freq)/num_chunks);
    freq_chunks = zeros([num_chunks, chunk_length]);
    p_freq = zeros([num_chunks, chunk_length]);
    r_squared_freq = zeros([num_chunks, 1]);
    freq_poly = zeros([num_chunks, 2]);

    for ii = 1:num_chunks
        freq_chunks(ii, :) = iq_freq((ii-1)*chunk_length+1 : ii*chunk_length);

        x_freq = (0:chunk_length-1) * (1/Fs);
        [freq_poly(ii, :), S_freq] = polyfit(x_freq, freq_chunks(ii, :), 1);
        ss_res = S_freq.normr ^ 2;
        ss_tot = sum((freq_chunks(ii, :) - mean(freq_chunks(ii, :))).^2);
        r_squared_freq(ii) = 1 - ss_res/ss_tot;
        p_freq(ii, :) = freq_poly(ii, 1) * x_freq + freq_poly(ii, 2);        
    end
    
    figure(200);
    plot(freq_chunks.');
    hold on;
    plot(p_freq(r_squared_freq > r_squared_threshold, :).', '.-');
    hold off;

    mean_slope = mean(freq_poly(r_squared_freq > r_squared_threshold, 1));
    if (isnan(mean_slope))
        %mean_slope = -20 * 1e3/1e-6; %500 Hz
        mean_slope = -8 * 1e3/1e-6; %50 Hz
    end

    t_freq = (0:length(iq_freq)-1).' * (1/Fs);
    
    f_chirp = -0.5 * mean_slope * t_freq;
    y_chirp = exp(1j*2*pi*f_chirp.*t_freq);
    
    iq_dechirped = iq_data(1:end-1) .* y_chirp;

    phase_dechirped = unwrap(atan2(imag(iq_dechirped), real(iq_dechirped)));
    f_dechirped = (1/(2*pi)) * diff(phase_dechirped) / (1/Fs);        
    

    figure(201);
    plot(t_freq, iq_freq, t_freq(1:end-1), f_dechirped);
    return;
end