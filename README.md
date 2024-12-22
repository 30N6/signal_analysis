# ARSR-4: Analysis of PDWs collected with pluto_esm
<img src="https://github.com/user-attachments/assets/87500873-1f59-4094-949e-ef7e03e506f0" width=20%>
<img src="https://github.com/user-attachments/assets/b6028ff0-afa0-477d-b4ee-cabb6d7d1f61" width=20%>
<img src="https://github.com/user-attachments/assets/43c5481e-17f0-4b02-97e4-dc6e5a06f665" width=20%>
<img src="https://github.com/user-attachments/assets/e69f4263-4187-4631-bfbe-59684714dd6d" width=20%>

## Introduction
This report presents a brief analysis of the ARSR-4 radar, based on PDWs collected with the pluto_esm open-source ESM receiver. Basic radar interpulse and intrapulse parameters are examined, including the amplitude/scan pattern, pulse duration, PRI, and modulation.

## Pluto_esm
The pluto_esm platform is an open-source implementation of a narrowband, channelized ESM receiver based on the ADALM-PLUTO evaluation board combined with Python software. 

On the hardware side, the FPGA is loaded with a modified ADALM-PLUTO image, where custom ESM logic has been added. This includes a line-rate 64-channel channelizer, a multi-channel PDW encoder, a dwell controller, reporting logic, etc.

The pluto_esm_app GUI controls the operation of the hardware, programming the FPGA with a frequency hopping program according to the user-specified scan configuration. For each hop/dwell, the FPGA reports channel statistics (accumulated and maximum power) for the 64 channelizer channels, and if desired, PDWs for any detected pulses. All reports are sent via DMA to the Linux OS running on the Zynq PS layer, where a small C program fetches them from memory and forwards the data to the USB-connected Ethernet adapter.

<img src="https://github.com/user-attachments/assets/87617445-a41e-4e44-8132-8c5be44d5f16" width=75%>

### Pluto_esm specifications
| Parameter            | Value                      |
|----------------------|:-----------:               |
| Transceiver          | AD9363                     |
| Tuning range         | 70-6000 MHz                | 
| Instantaneous BW     | 56 MHz                     |
| FPGA                 | XC7Z010                    |
| Data link            | 100M Ethernet (emulated)   |
| Channelizer channel count   | 64                  |
| Channelizer spacing (~BW)   | 0.96 MHz            |
| Channelizer sampling freq   | 1.92 MHz            |

### Pluto_esm PDWs
PDWs are stored by pluto_esm_app in json format. An example is annotated below:
```{"time": [2024, 12, 18, 15, 58, 29, 2, 353, 0], "sec_frac": 0.6321597099304199,   # logging timestamp
  "data": {
    "msg_seq_num": 0,               # internal message sequence number
    "msg_type": 32,                 # internal message type
    "dwell_seq_num": 0,             # dwell sequence number
    "pulse_seq_num": 7,             # per-channel pulse sequence number
    "pulse_channel": 37,            # channelizer channel index
    "pulse_threshold": 8,           # current threshold
    "pulse_power_accum": 1292,      # pulse power (I^2 + Q^2) accumulator
    "pulse_duration": 69,           # pulse duration, in channelizer channel cycles (1.92 MHz)
    "pulse_frequency": 0,           # currently unused
    "pulse_start_time": 320421454,  # pulse TOA, in system clock cycles (245.76 MHz, 4x the ADC sampling frequency)
    "buffered_frame_index": 0,      # internal IQ capture index
    "buffered_frame_valid": 1,      # IQ capture valid flag
    "buffered_frame_data": [[1, -1], [0, -1], ... ],  # raw IQ data of the pulse: 8 samples before the trigger point, 40 samples after
    "channel_frequency": 1252.8,    # channelizer channel frequency
    "dwell_channel_entry": {        # channel statistics (spectrum analyzer)
      "index": 37,                    # channelizer channel index
      "accum": 236740,                # channel power accumulator (I^2 + Q^2)
      "max": 102                      # channel power max value for the current dwell
    },
    "dwell_threshold_shift": 3,     # automatic threshold control setting
    "modulation_data": {                      # modulation analysis performed by pluto_esm_app using the raw IQ data
      "modulation_type": "FM",                  # frequency modulation - only LFM currently supported
      "LFM_slope": -15918.92858875841,          # calculated LFM slope, Hz/us
      "LFM_r_squared": 0.634858617417082,       # R^2 - goodness of fit of the calculated slope
      "LFM_mean_residual": 54475.58927354079    # mean residual of the slope fit
    }
  }
}
```

## Collection setup
* ARSR-4 to receiver: R=17 mi, $$\epsilon$$=-2.5&deg;
* Data collected using a vertically-polarized broadband (700-6000 MHz) antenna
* To improve coverage, an emitter-specific config is used, which specifies a 314 ms dwell (not divisible by the PRI or scan time) near each of the two radar frequencies
```
{
    "sim_mode": {"enable": 0, "filename": ""},
    "enable_recording": 0,
    "analysis_config": {
      "enable_pdw_recording": 1,
      "modulation_threshold": 0.25,
      "pulsed_emitter_search": {"expected_pulse_count": 0.1, "PW_range_scaling": [0.25, 1.25], "PRI_range_scaling": [0.75, 1.25]},
      "modulation_analysis": {"FM_threshold_residual": 0.05, "FM_threshold_r_squared": 0.5, "FM_threshold_slope": 1000, "FM_min_samples": 8}
    },
    "fast_lock_config": {"recalibration_interval": 600.0, "recalibration_pause": 2.0},
    "dwell_config": {"freq_start": 96.0, "freq_step": 48.0, "channel_step": 0.96},
    "scan_config": {
      "randomize_scan_order": 0,
      "include_freqs": [
          {"freq_range": [1250, 1255], "dwell_time": 0.31415926, "comment": "ARSR-4"},
          {"freq_range": [1334, 1338], "dwell_time": 0.31415926, "comment": "ARSR-4"}
      ],
      "exclude_freqs": []
    },
    "emitter_config": {
        "pulsed_emitters": [
            {"name": "ARSR-4",  "freq_range": [1200, 1400], "PW_range": [60, 90],     "PRI_range": [1500, 14000], "priority": 2, "threshold_dB": 9}
        ],
        "cw_emitters": []
    },
    "pluto_dma_reader_path": "../pluto_dma_reader/pluto_dma_reader",
    "pluto_credentials": {"username": "root", "password": "analog"},
    "graphics": {"fullscreen": 0, "noframe": 0}
}
```
## Collection analysis

### Frequency
* Pulses were detected in three channels, centered at 1252.80, 1253.76, and 1336.32 MHz.
* Pulses are detected in the adjacent channels at 1252.80 and 1253.76 MHz due to FMOP (the frequency moving from one channel to another within the pulse).
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_1.png)

### Pulse duration
* The upper frequency exhibits a narrow peak near 59 us, while the lower two frequencies have broad peaks centered near 70 and 38 us.
* The lower channels, being adjacent, have overlapping frequency frequency responses. Therefore, the total duration of the lower frequency pulse should be lower than 70+38 = 108 us.
* Overall, the collected pulse durations are consistent with the values published in open literature: 60 us for the high frequency and 90 us for the low frequency.
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_2.png)

### Scan
* A 12 second scan period is clearly apparent.
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_3.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_3_detail.png)

### PRI
* The ARSR-4 uses a PRI stagger.
* Computing a full PRI histogram (a histogram of the TOA differences between a given pulse and multiple subsequent pulses), we find a single prominent peak at the common stagger sum, around 41667 us.
* As expected, the PRI pattern is the same between frequencies.
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_4.png)
* With a first-level PRI histogram (TOA differences of adjacent pulses only), the stagger pattern is easier to see.
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_5.png)
* Finally, we can identify the PRIs via automatic clustering.
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_9.png)
* Applying a threshold to eliminate spurious values, we find that there are nine PRIs:
```
PRI clustering, freq=1336.32:
   1: median: 1957.29  N: 3272
   3: median: 3336.98  N:13396
   5: median: 4203.65  N: 3628
   6: median: 4901.56  N: 3557
   8: median: 2382.29  N: 3206
   9: median: 3604.17  N: 3577
  11: median: 2839.58  N: 3240
  12: median: 4540.10  N: 3562
  13: median: 3891.15  N: 3605
```

### Raster/PRI stagger
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_10.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_11.png)
```
Stagger pattern from raster, freq=1336.32:
   1: PRI: 3337.03 us   PRF: 299.67 Hz
   2: PRI: 2382.05 us   PRF: 419.81 Hz
   3: PRI: 3891.40 us   PRF: 256.98 Hz
   4: PRI: 4203.45 us   PRF: 237.90 Hz
   5: PRI: 3336.61 us   PRF: 299.71 Hz
   6: PRI: 2839.61 us   PRF: 352.16 Hz
   7: PRI: 3336.90 us   PRF: 299.68 Hz
   8: PRI: 4901.66 us   PRF: 204.01 Hz
   9: PRI: 3337.08 us   PRF: 299.66 Hz
  10: PRI: 1956.84 us   PRF: 511.03 Hz
  11: PRI: 3604.05 us   PRF: 277.47 Hz
  12: PRI: 4540.30 us   PRF: 220.25 Hz
```

### Modulation
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_21.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_22.png)

```
>> modulation_analysis
SNR > 30.0 dB:
example [1]: normr= 2520.4  r_squared=0.983  snr= 43.4  44.2 -- slope: -20161.6
example [2]: normr= 3935.3  r_squared=0.957  snr= 58.5  77.3 -- slope: -19892.0
example [3]: normr= 3070.1  r_squared=0.974  snr= 88.3  61.2 -- slope: -20111.7
example [4]: normr= 3922.7  r_squared=0.957  snr=114.5 130.0 -- slope: -19827.4
Hardware detection summary: num_pulses=1920  num_detected_FM=1920 (100.0%)  mean_r_squared=0.972 mean_slope=-20307.1 mean_residual=10699.4

SNR < 15.0 dB:
example [1]: normr= 6955.5  r_squared=0.860  snr=  3.2   3.4 -- slope: -18337.1
example [2]: normr= 8001.2  r_squared=0.802  snr=  2.4   1.4 -- slope: -17179.3
example [3]: normr= 8647.9  r_squared=0.841  snr=  2.6   8.0 -- slope: -21209.2
example [4]: normr=12752.5  r_squared=0.707  snr=  1.5   4.2 -- slope: -21109.1
Hardware detection summary: num_pulses=4106  num_detected_FM=4105 (100.0%)  mean_r_squared=0.806 mean_slope=-19975.4 mean_residual=44907.4
```
