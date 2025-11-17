import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

calibrations_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\calibration_crystal_0_new.h5"
scatter_folderpathA = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_scatters"
scatter_folderpathB = r"H:\alpha_spect_mini\20250225_pc123456_usb123456_daq123456_800fps_30min_Co57\panel_1\crystals_scatters"

with h5py.File(calibrations_filepath, 'r') as f:
    slopes = f['pixels_calibrations'][:, 0]
    intercepts = f['pixels_calibrations'][:, 1]
    isvalids = f['pixels_isvalids'][:]
    thresholds = f['pixels_thresholds'][:]
    peaks = f['pixels_peaks'][:]



for pixel_id in range(0, 6400):
    row = pixel_id // 80
    col = pixel_id % 80
    
    scatter_filepath = scatter_folderpathA + rf"\scatter_crystal_0_pixel_{pixel_id}.bin"
    if not Path(scatter_filepath).exists():
        continue
    scattersA = np.fromfile(scatter_filepath, dtype=np.uint16)

    scatter_filepath = scatter_folderpathB + rf"\scatter_crystal_0_pixel_{pixel_id}.bin"
    if not Path(scatter_filepath).exists():
        continue
    scattersB = np.fromfile(scatter_filepath, dtype=np.uint16)

    plt.figure(figsize=(10, 6))
    plt.subplot(1, 2, 1)
    plt.plot(scattersA, color='red', linestyle='none', marker='.', markersize=1)
    plt.plot(scattersB, color='blue', linestyle='none', marker='.', markersize=1)
    plt.axhline(thresholds[pixel_id], color='blue', linestyle='--')
    plt.axhline(peaks[pixel_id, 0], color='green', linestyle='--')
    plt.ylabel("ADU")

    plt.subplot(1, 2, 2)
    erg_scattersA = scattersA * slopes[pixel_id] + intercepts[pixel_id]
    erg_scattersB = scattersB * slopes[pixel_id] + intercepts[pixel_id]
    plt.plot(erg_scattersA, color='red', linestyle='none', marker='.', markersize=1)
    plt.plot(erg_scattersB, color='blue', linestyle='none', marker='.', markersize=1)
    plt.axhline(thresholds[pixel_id] * slopes[pixel_id] + intercepts[pixel_id], color='blue', linestyle='--')
    plt.axhline(peaks[pixel_id, 0] * slopes[pixel_id] + intercepts[pixel_id], color='green', linestyle='--')
    plt.ylabel("Energy (keV)")

    plt.suptitle("Pixel Thresholds on Scatter Plot - Panel 1, Crystal 0, Pixel {}".format(pixel_id))
    plt.show()