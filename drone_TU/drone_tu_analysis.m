file_dir = "./data";
drone_entries = {{'drone_tu_preprocessed_DJI_inspire_2_2G.mat',               2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_inspire_2_5G_1of2.mat',          5.800e9, 200e6};
                 {'drone_tu_preprocessed_DJI_inspire_2_5G_2of2.mat',          5.800e9, 200e6};
                 {'drone_tu_preprocessed_DJI_matrice_100_2G.mat',             2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_matrice_210_2G.mat',             2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_matrice_210_5G_2of2.mat',        5.800e9, 200e6};
                 {'drone_tu_preprocessed_DJI_mavic_mini_2G.mat',              2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_mavic_pro_2G.mat',               2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_phantom_4_2G.mat',               2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_phantom_4_pro_plus_2G.mat',      2.440e9, 120e6};
                 {'drone_tu_preprocessed_DJI_phantom_4_pro_plus_5G_1of2.mat', 5.800e9, 200e6};
                 {'drone_tu_preprocessed_DJI_phantom_4_pro_plus_5G_2of2.mat', 5.800e9, 200e6};
                 {'drone_tu_preprocessed_Parrot_disco_2G.mat',                2.440e9, 120e6};
                 {'drone_tu_preprocessed_Parrot_mambo_control_2G.mat',        2.440e9, 120e6};
                 {'drone_tu_preprocessed_Parrot_mambo_video_2G.mat',          2.440e9, 120e6};
                 {'drone_tu_preprocessed_Yuneec_typhoon_h_2G_1of2.mat',       2.440e9, 120e6};
                 {'drone_tu_preprocessed_Yuneec_typhoon_h_2G_2of2.mat',       2.440e9, 120e6};
                 {'drone_tu_preprocessed_Yuneec_typhoon_h_5G.mat',            5.700e9, 200e6};
                };

%for ii = 1:length(drone_entries)
ii = 18;

filename = file_dir + "/" + drone_entries{ii}{1};
d = load(filename);

dt = 1/d.input_fs;

padded_len = 2^ceil(log2(length(d.iq_data_m)));
iq_data = paddata(d.iq_data_m, padded_len);
iq_phase = unwrap(atan2(imag(iq_data), real(iq_data)));
iq_freq = (1/(2*pi)) * diff(iq_phase) / dt; 

reports = struct("iq_data", [], "iq_phase", [], "iq_freq", [], "iq_length", 0);

for ii = 1:size(iq_data, 2)
    reports(ii).fs = d.input_fs;
    reports(ii).iq_data = iq_data(:, ii);
    reports(ii).iq_phase = iq_phase(:, ii);
    reports(ii).iq_freq = [iq_freq(:, ii); iq_freq(end, ii)];
    reports(ii).iq_length = size(iq_data, 1);
end



%%
offset = 0;
num_rows = 4;
num_cols = 4;

figure(2);
ax1 = zeros(num_rows, num_cols);
figure(3);
ax2 = zeros(num_rows, num_cols);

for row = 1:num_rows
    for col = 1:num_cols
        plot_index = (row-1) * num_cols + col;
        d = reports(plot_index + offset);        
        
        t = (0:d.iq_length-1).' * dt;
        y = d.iq_data;
        
        figure(2);
        ax1(row, col) = subplot(num_rows,num_cols,plot_index);
        plot(t, real(y), t, imag(y));
    end
end

for row = 1:num_rows
    for col = 1:num_cols
        plot_index = (row-1) * num_cols + col;
        d = reports(plot_index + offset);        

        t = (0:d.iq_length-1).' * dt;
        y = d.iq_freq;

        mean_freq = mean(y);
        
        figure(3);
        ax1(row, col) = subplot(num_rows,num_cols,plot_index);
        plot(t, y, [t(1), t(end)], [mean_freq, mean_freq]);
    end
end
%%

d = reports(3);

t = (0:length(d.iq_data)-1).' * dt;

% d.iq_data = d.iq_data(t < 0.5e-4);
% t = t(t < 0.5e-4);

[fft_median, fft_mean, fft_std] = get_fft_stats(d.iq_data, d.fs)

%mean_freq = mean(d.iq_freq);
%mean_freq = fft_mean;
mean_freq = fft_median;

shift_freq = -mean_freq - 0; %.25e6;


iq_shifted_data = d.iq_data .* exp(2j*pi*shift_freq.*t);
iq_shifted_phase = unwrap(atan2(imag(iq_shifted_data), real(iq_shifted_data)));
iq_shifted_freq = (1/(2*pi)) * diff(iq_shifted_phase) / dt; 

lpf = fir1(1024, pi*2e6/d.fs);

iq_filtered_data = filter(lpf, 1, iq_shifted_data);
iq_filtered_phase = unwrap(atan2(imag(iq_filtered_data), real(iq_filtered_data)));
iq_filtered_freq = (1/(2*pi)) * diff(iq_filtered_phase) / dt; 
%iq_shifted_phase = unwrap(atan2(imag(iq_shifted_data), real(iq_shifted_data)));
%iq_shifted_freq = (1/(2*pi)) * diff(iq_shifted_phase) / dt; 

fs_d = d.fs/16;
iq_ds_data = downsample(iq_filtered_data,16);
t_ds = (0:length(iq_ds_data)-1).' * (1/fs_d);
%t = (0:d.iq_length-1).' * dt;
%iq_downsampled = decimate(iq_shifted_data, 16);


L = length(d.iq_data);
f_shifted = (d.fs/L)*(-L/2:L/2-1).';

Y_orig = fft(d.iq_data);
Y_orig_shifted = fftshift(Y_orig);
Y_orig_abs = abs(Y_orig_shifted);

Y_shifted = fft(iq_shifted_data);
Y_shifted_shifted = fftshift(Y_shifted);
Y_shifted_abs = abs(Y_shifted_shifted);

Y_filtered = fft(iq_filtered_data);
Y_filtered_shifted = fftshift(Y_filtered);
Y_filtered_abs = abs(Y_filtered_shifted);


figure(5);
subplot(4,1,1);
plot(t(1:end-1), iq_shifted_freq, t(1:end-1), iq_filtered_freq);

subplot(4,1,2);
[xc, lags] = xcorr(d.iq_data);
plot(lags*dt, abs(xc));

subplot(4,1,3);
[xc, lags] = xcorr(iq_ds_data);
plot(lags*(1/fs_d), abs(xc));

subplot(4,1,4);
N = 32;
r = compute_sfft(iq_ds_data, N);
imagesc(r.');


% figure(6);
% subplot(4,1,1);
% instfreq(d.iq_data, d.fs);
% subplot(4,1,2);
% instfreq(iq_shifted_data, d.fs);
% subplot(4,1,3);
% instfreq(iq_filtered_data, d.fs);
% subplot(4,1,4);
% instfreq(iq_ds_data, fs_d);

figure(7);
subplot(4,1,1);
plot(t, iq_shifted_data, t, iq_filtered_data, t_ds, iq_ds_data, '.-');

subplot(4,1,2);
max_Y = max(Y_orig_abs);
plot(f_shifted, Y_orig_abs, f_shifted, Y_shifted_abs, f_shifted, Y_filtered_abs, [fft_mean, fft_mean], [0, max_Y]);

subplot(4,1,3);
plot(t, iq_shifted_phase);
subplot(4,1,4);
plot(t(1:end-1), iq_shifted_freq);


function r = compute_sfft(data, N)
    padded_data = paddata(data, 2^ceil(log2(length(data))));

    r = zeros(length(padded_data)/N, N);
    for ii = 1:(length(padded_data)/N)
        r(ii, :) = fftshift(abs(fft(padded_data((ii-1)*N+1 : ii*N))));
    end
end


function [freq_median, freq_mean, freq_std] = get_fft_stats(data, Fs)
    L = length(data);
    f_shifted = (Fs/L)*(-L/2:L/2-1).';

    Y = fft(data);
    Y_shifted = fftshift(Y);
    
    Y_abs = abs(Y_shifted);
    Y_sum = sum(Y_abs);
    Y_scaled = Y_abs .* (1/Y_sum);
    Y_weighted_1 = Y_scaled .* f_shifted;
    
    Y_cs = cumsum(Y_scaled);
    i_first = find(Y_cs >= 0.5, 1);
    freq_median = f_shifted(i_first);

    freq_mean = sum(Y_weighted_1);
    Y_var = sum(Y_scaled .* (f_shifted - freq_mean).^2);

    freq_std = sqrt(Y_var);
end


