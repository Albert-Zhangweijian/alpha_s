import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

calibrations_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\calibration_crystal_0_new.h5"
scatter_folderpath = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_scatters"
scatter_folderpath = r"H:\alpha_spect_mini\20250225_pc123456_usb123456_daq123456_800fps_30min_Co57\panel_1\crystals_scatters"

with h5py.File(calibrations_filepath, 'r') as f:
    slopes = f['pixels_calibrations'][:, 0]
    intercepts = f['pixels_calibrations'][:, 1]
    isvalids = f['pixels_isvalids'][:]
    thresholds = f['pixels_thresholds'][:]
    peaks = f['pixels_peaks'][:]

avg_adus = np.zeros(6400)
for pixel_id in range(6400):
    scatter_filepath = scatter_folderpath + rf"\scatter_crystal_0_pixel_{pixel_id}.bin"
    if not Path(scatter_filepath).exists():
        continue
    scatters = np.fromfile(scatter_filepath, dtype=np.uint16)
    avg_adus[pixel_id] = np.mean(scatters)
nondata_pixels = avg_adus == 0

avg_adu_minus_threshold = avg_adus - thresholds
avg_adu_minus_threshold[nondata_pixels] = 0

plt.imshow(avg_adu_minus_threshold.reshape(80, 80), cmap='gray')
plt.title("Average ADU Scatter Plot - Panel 1, Crystal 0")
plt.colorbar(label='ADU - Threshold ADU')
plt.clim(-100, 100)
plt.show()


for pixel_id in range(0, 6400):
    row = pixel_id // 80
    col = pixel_id % 80
    
    scatter_filepath = scatter_folderpath + rf"\scatter_crystal_0_pixel_{pixel_id}.bin"
    if not Path(scatter_filepath).exists():
        continue
    scatters = np.fromfile(scatter_filepath, dtype=np.uint16)

    plt.figure(figsize=(10, 6))
    plt.subplot(1, 2, 1)
    plt.plot(scatters, color='red', linestyle='none', marker='.', markersize=1)
    plt.axhline(thresholds[pixel_id], color='blue', linestyle='--')
    plt.axhline(peaks[pixel_id, 0], color='green', linestyle='--')
    plt.title("Pixel Thresholds on Scatter Plot - Panel 1, Crystal 0")
    plt.xlabel("Pixel Column")
    plt.ylabel("Pixel Row")

    plt.subplot(1, 2, 2)
    erg_scatters = scatters * slopes[pixel_id] + intercepts[pixel_id]
    plt.plot(erg_scatters, color='red', linestyle='none', marker='.', markersize=1)
    plt.axhline(thresholds[pixel_id] * slopes[pixel_id] + intercepts[pixel_id], color='blue', linestyle='--')
    plt.axhline(peaks[pixel_id, 0] * slopes[pixel_id] + intercepts[pixel_id], color='green', linestyle='--')
    plt.title("Pixel Thresholds on Scatter Plot (Energy) - Panel 1, Crystal 0")
    plt.xlabel("Pixel Column")
    plt.ylabel("Energy (keV)")

    plt.show()