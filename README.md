

![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_1.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_2.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_3.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_3_detail.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_4.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_5.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_9.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_10.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_11.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_21.png)
![image](https://github.com/30N6/radar_analysis/blob/master/ARSR_4/analysis-20241218-155827-ARSR-4_fig_22.png)

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
