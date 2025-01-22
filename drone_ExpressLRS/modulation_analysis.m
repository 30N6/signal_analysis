channel_clock_period    = (32/61.44e6);

pdw_reports = pdw_reports(([pdw_reports.buffered_frame_valid] == 1));

%plot([matching_pdws.pulse_toa_s], 10*log10([matching_pdws.implied_pulse_snr]), [matching_pdws.pulse_toa_s], 10*log10([matching_pdws.recorded_pulse_snr]), [matching_pdws.pulse_toa_s], 10*log10([matching_pdws.pulse_power]));

pdw_desc = {"SNR > 20.0 dB"; "SNR < 15.0 dB"};
pdws_by_snr = {pdw_reports(10*log10([pdw_reports.pulse_power]) > 30.0); pdw_reports(10*log10([pdw_reports.pulse_power]) < 15.0)};

num_subplots = 12;
dt = channel_clock_period;

for jj = 1:length(pdws_by_snr)
    figure(20 + jj);
    matching_pdws = pdws_by_snr{jj};

    if length(matching_pdws) < 1
        continue
    end
    
    fprintf('%s:\n', pdw_desc{jj});
    for ii = 1:num_subplots
        pdw = matching_pdws(ii);
        
        last_iq_index = min(8 + pdw.pulse_duration, 48);

        phase = unwrap(atan2(imag(pdw.recorded_iq_data(9:last_iq_index)), real(pdw.recorded_iq_data(9:last_iq_index))));
        freq = (1/(2*pi)) * diff(phase) / dt;
        
        x_phase = (0:length(phase)-1) * dt;
        [p_p, S_phase] = polyfit(x_phase, phase, 2);
        ss_res = S_phase.normr ^ 2;
        ss_tot = sum((freq - mean(freq)).^2);
        r_squared_phase = 1 - ss_res/ss_tot;
        p_phase = p_p(1) * (x_phase.^2) + p_p(2) * x_phase + p_p(3);
    
        x_freq = (0:length(freq)-1) * dt;
        [p_f, S_freq] = polyfit(x_freq, freq, 1);
        ss_res = S_freq.normr ^ 2;
        ss_tot = sum((freq - mean(freq)).^2);
        r_squared_freq = 1 - ss_res/ss_tot;
        p_freq = p_f(1) * x_freq + p_f(2);
        
        fprintf("example [%d]: normr=%7.1f  r_squared=%.3f  snr=%5.1f %5.1f -- slope: %.1f\n", ii, S_freq.normr / length(x_freq), r_squared_freq, ...
            pdw.implied_pulse_snr, pdw.recorded_pulse_snr, p_f(1) * 1e-6);
    
        subplot(num_subplots, 3, 3 * ii - 2);
        plot(1:length(pdw.recorded_iq_data(1:last_iq_index)), real(pdw.recorded_iq_data(1:last_iq_index)), ...
             1:length(pdw.recorded_iq_data(1:last_iq_index)), imag(pdw.recorded_iq_data(1:last_iq_index)), ...
             1:length(pdw.recorded_iq_data(1:last_iq_index)), abs(pdw.recorded_iq_data(1:last_iq_index)));
        grid on;
        if (ii == 1)
            title('IQ, amplitude');
        end
        ylabel('amplitude');
        if (ii == num_subplots)
            xlabel('Sample index');
        end
        
        subplot(num_subplots, 3, 3 * ii - 1);
        plot(x_phase, phase); %, x_phase, p_phase);
        grid on;
        if (ii == 1)
            title('Unwrapped phase');
        end
        ylabel('Phase (rad)');
        if (ii == num_subplots)
            xlabel('Time (s)');
        end
    
        subplot(num_subplots, 3, 3 * ii - 0);
        plot(x_freq, freq, 'o-', x_freq, p_freq);
        grid on;
        if (ii == 1)
            title('Frequency');
        end
        ylabel('Frequency (Hz)');
        if (ii == num_subplots)
            xlabel('Time (s)');
        end    
    end
    
    num_pulses      = length(matching_pdws);
    num_detected_FM = 0;
    sum_r_squared   = 0;
    sum_slope       = 0;
    sum_residual    = 0;
    for ii = 1:num_pulses
        pdw = matching_pdws(ii);
        if ~isstruct(pdw.modulation_data)
            continue;
        end

        %disp(pdw.modulation_data)
        if isstring(pdw.modulation_data.modulation_type) && strcmp(pdw.modulation_data.modulation_type, 'FM')
            num_detected_FM = num_detected_FM + 1;
            sum_r_squared   = sum_r_squared + pdw.modulation_data.LFM_r_squared;
            sum_slope       = sum_slope + pdw.modulation_data.LFM_slope * -1; %2024-12-29: old version had I and Q swapped
            sum_residual    = sum_residual + pdw.modulation_data.LFM_mean_residual;
        end
    end
    fprintf('Hardware detection summary: num_pulses=%d  num_detected_FM=%d (%.1f%%)  mean_r_squared=%.3f mean_slope=%.1f mean_residual=%.1f\n\n', ...
        num_pulses, num_detected_FM, 100 * (num_detected_FM/num_pulses), sum_r_squared/num_detected_FM, sum_slope/num_detected_FM, sum_residual/num_detected_FM);
    
end

figs_to_save = [21, 22];
split_fn = split(filename, '.');
filename_base = split_fn{1};

for ii = 1:length(figs_to_save)
    fig_index = figs_to_save(ii);
    f = figure(fig_index);
    f.Position = [500 100 800 1200];
    
    fig_filename = sprintf('%s_fig_%d.png', filename_base, fig_index);
    saveas(f, fig_filename);
end