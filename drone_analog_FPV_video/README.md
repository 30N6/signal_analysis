# Analog FPV drone video (NTSC): pluto_esm collection and analysis
<img src="https://github.com/user-attachments/assets/71724dc5-f3eb-4f7b-ba00-46dd7bb825bc" width=20%>

<img src="https://github.com/user-attachments/assets/6173134b-550b-4f20-bed6-bb3213df2b7f" width=20%>
<img src="https://github.com/user-attachments/assets/0ff5af5a-fadf-48ee-a5e9-540346ccfb84" width=20%>
<img src="https://github.com/user-attachments/assets/5ea475b8-3a45-4a3c-bedc-c3a291278b7d" width=20%>

## Introduction
This report presents a brief look at the NTSC video signal, commonly used on FPV drones, collected with the pluto_esm system (https://github.com/30N6/sw/wiki/Pluto_esm).

Although the NTSC waveform is complex and not particularly well-suited for a pulse-domain analysis (unlike, say, a radar), it does feature some periodic AM characteristics which allow a pulse-oriented ESM receiver to receive it as such. In particular, the fixed 29.97 fps frame rate together with the repeating 525 scan line/frame pattern provides a relatively stable PRI pattern among the detected pulses, and pluto_esm can effectively discriminate it in a noisy environment even with short dwells.

<img src="https://github.com/user-attachments/assets/9c9ea545-7fc2-4c0f-9644-6d51512e9c51" width=50%>


## Setup
* An inexpensive 5.8 GHZ FPV camera was used to generate the signal, with the power set to 20 mW.
* A number of pluto_esm configurations were tried, and it was found that even a short 10 ms dwell time setting performed well.
* Collection was performed indoors, with the antenna was placed 40 feet away from the transmitter in another room, in a somewhat congested RF environment.

<img src="https://github.com/user-attachments/assets/ec2f298d-cba2-460a-bca2-253baa05be3f" width=20%>


## Pluto_esm results
[PDW data](./analysis-20250118-143522-5917-20mW-10ms.log)

Starting with the full 5.6-5.95 GHz data, we can look at the pulse count by frequency. 
![image](./analysis-20250118-143522-5917-20mW-10ms_fig_1.png)

There is a large peak around 5917 MHz, the signal of interest. Large numbers of pulses from other sources appear across the band, but these aren't interesting and are excluded from further analysis (they do not preclude pluto_esm from detecting and identifying video signals). Narrowing the dataset to frequencies between 5900 and 5950 MHz, and excluding channels with fewer than 2000 pulses, there are three frequencies remaining in this collection. 

The pulse durations are mostly very short, 1-3 IQ samples long.
![image](./analysis-20250118-143522-5917-20mW-10ms_fig_2.png)

Pulse power varies randomly without a discernable pattern. Note that the peak to average power ratio is generally 3-6 dB, which limits detection range.
![image](./analysis-20250118-143522-5917-20mW-10ms_fig_3.png)

![image](./analysis-20250118-143522-5917-20mW-10ms_fig_5.png)
![image](./analysis-20250118-143522-5917-20mW-10ms_fig_10.png)
![image](./analysis-20250118-143522-5917-20mW-10ms_fig_11.png)

## Spectrum analyzer comparison
* Zero span mode not helpful - swamped with WiFi
![Screenshot 2025-01-18 143408- FPV video 5917 - PRI-based detection](https://github.com/user-attachments/assets/06701dd0-be4f-48aa-a597-b2b96caf22f2)

![image](https://github.com/user-attachments/assets/7978fc9a-2960-425e-994b-e9832fb3aa7f)

![image](https://github.com/user-attachments/assets/8b0eff8c-d231-4ca0-8127-ed8c4f545553)
![image](https://github.com/user-attachments/assets/078538c9-ccef-49a1-abf2-902c89ebc9c4)
![image](https://github.com/user-attachments/assets/008f8ad4-5d5c-410d-ba81-6c73e486cdee)



## Conclusion
Although pluto_esm relies on peak power for detection (with no integration across pulses), its sensitivity appears to be better than that of a similarly-configured spectrum analyzer with this type of signal. 


