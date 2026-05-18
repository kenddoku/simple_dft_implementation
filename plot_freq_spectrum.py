import numpy as np
import matplotlib.pyplot as plt
import math

with open("signal.txt", "r") as f:
    f.readline()
    f.readline()

    signal_data = np.loadtxt(f)
    signal_time, signal = signal_data[:,0], signal_data[:,1]

with open("samples.txt", "r") as f:
    # line = f.readline() # "sine_freqs[Hz]: <sine_freq0> <sine_freq1> [...]\n"
    # sine_freqs      = np.array((line.split(":")[1]).split(" ")[1:], dtype=float)
    
    # line = f.readline() # "sine_phases[rad]: <init_phase0> <init_phase1> [...]\n" 
    # initial_phase   = np.array((line.split(":")[1]).split(" ")[1:], dtype=float)
    
    # line = f.readline() # "sampling_freq[Hz]: <sampling_freq>\n"
    # sampling_freq   = float(line.split(":")[1])
    
    # line = f.readline() # "sampling_time[s]: <sampling_time>\n"
    # sampling_time   = float(line.split(":")[1])
    
    # line = f.readline() # "samples_len: <samples_len>\n"
    # samples_len     = int(line.split(":")[1])

    samples_data = np.loadtxt(f)
    samples_time, samples = samples_data[:,0], samples_data[:,1]

# print(f"sampling_freq[Hz]: {sampling_freq}")
# print(f"sampling_time[s]: {sampling_time}")
# print(f"samples_len: {samples_len}")

with open("freq_magnitudes_phases.txt") as file:
    freq_magnitudes_phases = np.loadtxt(file)
freq, freq_mag, freq_phase = freq_magnitudes_phases[:,0], freq_magnitudes_phases[:,1], freq_magnitudes_phases[:,2]

# Plotting results ==================================================================
fig, axes = plt.subplots(2, 1)
axes[0].plot(signal_time, signal, c='b', label="signal", zorder=1)
axes[0].scatter(samples_time, samples, s=1.5, c='r', label="samples", zorder=2)
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel("amplitude [V]") 
axes[0].set_title("Original signal (blue) and sampled points (red)")

axes[1].plot(freq, freq_mag)
axes[1].set_xlabel("frequency [Hz]")
axes[1].set_ylabel("freq_magnitude")
axes[1].set_title("Frequency spectrum")

plt.tight_layout()
plt.show()
