data_path = 'C:/drone_data/Radio-Frequency Control and Video Signal Recordings of Drones';

% The drones were recorded at 2.44 GHz and 5.8 GHz center frequencies (with the exception of Yuneec Typhoon H at 5.7 GHz) if the drone supported both.
% In the dataset, these center frequencies were denoted with 2G and 5G in the filename, respectively.
%The sampling frequency was 120 MHz at 2.44 GHz and 200 MHz at 5.8 GHz.


reload = 0;
N_filt = 1024;
t_segment = 250e-6;
num_segments = 64;

drone_entries = {{'DJI_inspire_2_2G.bin',               2.440e9, 120e6, 0.1};
                 {'DJI_inspire_2_5G_1of2.bin',          5.800e9, 200e6, 0.1};
                 {'DJI_inspire_2_5G_2of2.bin',          5.800e9, 200e6, 0.1};
                 {'DJI_matrice_100_2G.bin',             2.440e9, 120e6, 0.1};
                 {'DJI_matrice_210_2G.bin',             2.440e9, 120e6, 0.001};
                 {'DJI_matrice_210_5G_2of2.bin',        5.800e9, 200e6, 0.01};
                 {'DJI_mavic_mini_2G.bin',              2.440e9, 120e6, 0.1};
                 {'DJI_mavic_pro_2G.bin',               2.440e9, 120e6, 0.1};
                 {'DJI_phantom_4_2G.bin',               2.440e9, 120e6, 0.1};
                 {'DJI_phantom_4_pro_plus_2G.bin',      2.440e9, 120e6, 0.1};
                 {'DJI_phantom_4_pro_plus_5G_1of2.bin', 5.800e9, 200e6, 0.1};
                 {'DJI_phantom_4_pro_plus_5G_2of2.bin', 5.800e9, 200e6, 0.1};
                 {'Parrot_disco_2G.bin',                2.440e9, 120e6, 0.1};
                 {'Parrot_mambo_control_2G.bin',        2.440e9, 120e6, 0.1};
                 {'Parrot_mambo_video_2G.bin',          2.440e9, 120e6, 0.1};
                 {'Yuneec_typhoon_h_2G_1of2.bin',       2.440e9, 120e6, 0.1};
                 {'Yuneec_typhoon_h_2G_2of2.bin',       2.440e9, 120e6, 0.1};
                 {'Yuneec_typhoon_h_5G.bin',            5.700e9, 200e6, 0.1};
                };

for i_entry = 1:length(drone_entries)
    filename = drone_entries{i_entry}{1};
    input_f0 = drone_entries{i_entry}{2};
    input_fs = drone_entries{i_entry}{3};
    input_threshold = drone_entries{i_entry}{4};
    
    fprintf("Loading %s\n", filename);

    iq_data = get_drone_data(data_path + "/" + filename);
    dt = 1/input_fs;
    t = (0:length(iq_data)-1) * dt;    
    %iq_phase = unwrap(atan2(imag(iq_data), real(iq_data)));
    %iq_freq = (1/(2*pi)) * diff(iq_phase) / dt; 
    iq_power = abs(iq_data).^2;

    N_segment = t_segment * input_fs;
    
    iq_power_avg = filter(ones([N_filt,1])/N_filt, 1, iq_power);
    
    iq_power_thresh = iq_power_avg > input_threshold;
    d_pt = diff(iq_power_thresh) > 0;
    d_pt_i = find(d_pt);
    d_pt_i = d_pt_i(2:end) - N_filt/2;

    output_len = min(num_segments, length(d_pt_i) - 2);
    iq_data_m = zeros([N_segment, output_len]);
    
    for ii = 1:output_len
        i_segment = d_pt_i(ii);
        iq_data_m(:, ii) = iq_data(i_segment:(i_segment + N_segment - 1));
    end

    t_dpt = t(d_pt_i(1:output_len));
    y_dpt = ones([output_len, 1]);
    t_slice = t(1:d_pt_i(output_len));
    y_slice = iq_data(1:d_pt_i(output_len));
    p_slice = iq_power_avg(1:d_pt_i(output_len));
    figure(1);
    plot(t_slice, real(y_slice), t_slice, imag(y_slice), t_slice, p_slice, t_dpt, y_dpt, 'o');

    split_filename = split(filename, ".");
    output_filename = sprintf("./data/drone_tu_preprocessed_%s.mat", split_filename{1});

    fprintf("Saving %s\n", output_filename);

    save(output_filename, "iq_data_m", "filename", "input_f0", "input_fs");
end

function iq_data = get_drone_data(filepath)
    iq_data = single(drone_data_tu_load_bin(filepath));
    iq_data = (iq_data-mean(iq_data))/(sqrt(var(iq_data)));
end
