import h5py
import numpy as np
import matplotlib.pyplot as plt
from lib.interactive_map import *


"""
This script is to show the thresholds value calculated from percentiles, the unit is keV.
"""


panel_id = 1
crystal_id = 0
calibration_file = rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\calibration_crystal_{crystal_id}_new.h5"
threshold_file = rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\pixels_thresholds_crystal_{crystal_id}_95per.h5"



# load the calibration params and 5keV thresholds
with h5py.File(calibration_file, 'r') as f:
    pixels_calibrations = f['pixels_calibrations'][:]
    pixels_thresholds_default = f['pixels_thresholds'][:]
    pixels_isvalids = f['pixels_isvalids'][:]

# load the measured thresholds by percentiles
with h5py.File(threshold_file, 'r') as f:
    pixels_thresholds_mesured_ADU = f['pixels_thresholds'][:]
    pixels_thresholds_measured = pixels_thresholds_mesured_ADU * pixels_calibrations[:, 0] + pixels_calibrations[:, 1]



pixel1_threshold_default = pixels_calibrations[0,0]*pixels_thresholds_default[0] + pixels_calibrations[0,1]
print(f"For panel {panel_id}, crystal {crystal_id}: default threshold={pixel1_threshold_default:.2}keV, measured threshold={pixels_thresholds_measured[0]:.2}keV")


# filter bad pixels
pixels_thresholds_measured[~pixels_isvalids] = 0
pixels_thresholds_default[~pixels_isvalids] = 0


pixels_spectra_file = rf"H:\alpha_spect_mini\nosource_500fps_3hours\panel_{panel_id}\pixels_spectra_crystal_{crystal_id}.h5"
with h5py.File(pixels_spectra_file, 'r') as f:
    pixels_spectra = np.array(f['pixels_spectra'][:])


r = ThresholdMap(pixels_thresholds_measured, pixels_spectra, maprange=(0, 20), shape=(80, 80), 
             thresholds=pixels_thresholds_mesured_ADU, adu_bins=None, calibration_params=pixels_calibrations)
r.spectrum_range = (0, 20)
plt.title(f"Panel {panel_id}, Crystal {crystal_id} - Thresholds Map")
plt.show()