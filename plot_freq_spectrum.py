import numpy as np
import matplotlib.pyplot as plt
import math

with open("sampled_sine_wave.txt", "r") as f:
    line = f.readline()         # "sine_freq[Hz]: <sine_freq>\n"
    sine_freq       = float(line.split(":")[1])
    
    line = f.readline()
    initial_phase   = float(line.split(":")[1])
    
    line = f.readline()         # "sampling_freq[Hz]: <sampling_freq>\n"
    sampling_freq   = float(line.split(":")[1])
    
    line = f.readline()         # "sampling_time[s]: <sampling_time>\n"
    sampling_time   = float(line.split(":")[1])
    
    line = f.readline()         # "samples_len: <samples_len>\n"
    samples_len     = int(line.split(":")[1])

    samples_data    = np.loadtxt(f)
    sampled_time, samples = samples_data[:,0], samples_data[:,1]

freq_ratio = sine_freq / sampling_freq

print(f"sine_freq[Hz]: {sine_freq}")   
print(f"initial_phase[rad]: {initial_phase}")
print(f"sampling_freq[Hz]: {sampling_freq}")
print(f"sampling_time[s]: {sampling_time}")
print(f"samples_len: {samples_len}")

with open("freq_magnitudes_phases.txt") as file:
    freq_magnitudes_phases = np.loadtxt(file)
freq, freq_mag, freq_phase = freq_magnitudes_phases[:,0], freq_magnitudes_phases[:,1], freq_magnitudes_phases[:,2]

time_original_signal = np.linspace(0.0, sampling_time, 10*samples_len)

fig, axes = plt.subplots(2, 1)
axes[0].plot(time_original_signal, np.sin(2*math.pi * time_original_signal * sine_freq + initial_phase), c='b', label="signal", zorder=1)
axes[0].scatter(sampled_time, samples, s=2, c='r', label="samples", zorder=2)
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel("amplitude [V]") 
axes[0].set_title("Original signal (blue) and sampled points (red)")

axes[1].plot(freq, freq_mag)
axes[1].set_xlabel("frequency [Hz]")
axes[1].set_ylabel("freq_magnitude")
axes[1].set_title("Frequency spectrum")

plt.tight_layout()
plt.show()
