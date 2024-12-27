analysis_freq           = 2893.44;
channel_clock_period    = (32/61.44e6);

reports = pdw_reports_by_freq{analysis_freq};
reports = reports(([reports.buffered_frame_valid] == 1));

%plot([matching_pdws.pulse_toa_s], 10*log10([matching_pdws.implied_pulse_snr]), [matching_pdws.pulse_toa_s], 10*log10([matching_pdws.recorded_pulse_snr]), [matching_pdws.pulse_toa_s], 10*log10([matching_pdws.pulse_power]));

pdw_desc = {"SNR > 30.0 dB"; "SNR < 15.0 dB"};
pdws_by_snr = {reports(10*log10([reports.pulse_power]) > 30.0); reports(10*log10([reports.pulse_power]) < 15.0)};

num_subplots = 4;
dt = channel_clock_period;

for jj = 1:length(pdws_by_snr)
    figure(20 + jj);
    matching_pdws = pdws_by_snr{jj};
    
    fprintf('%s:\n', pdw_desc{jj});
    for ii = 1:num_subplots
        pdw = matching_pdws(ii);
      
        x_i = 1:(pdw.pulse_duration + 8);

        subplot(num_subplots, 1, ii);
        plot(x_i, real(pdw.recorded_iq_data(x_i)), x_i, imag(pdw.recorded_iq_data(x_i)), x_i, abs(pdw.recorded_iq_data(x_i)));
        grid on;
        if (ii == 1)
            title('IQ, amplitude');
        end
        ylabel('amplitude');
        if (ii == num_subplots)
            xlabel('Sample index');
        end
    end
end

figs_to_save = [21, 22];
split_fn = split(filename, '.');
filename_base = split_fn{1};

for ii = 1:length(figs_to_save)
    fig_index = figs_to_save(ii);
    f = figure(fig_index);
    f.Position = [500 100 800 800];
    
    fig_filename = sprintf('%s_fig_%d.png', filename_base, fig_index);
    saveas(f, fig_filename);
end