import pathlib
import h5py
import numpy as np
import matplotlib.pyplot as plt
from sympy import false

from pyscripts.lib.calculate_thresholds import calculate_thresholds
from pyscripts.lib.cmd_wrapper import *
from pyscripts.lib.find_peaks import find_peaks
from pyscripts.lib.fit_peaks import fit_peaks, asymmetric_gaussian
from pyscripts.lib.calibrate_peaks import calibrate_peaks

#%% raw file processing

# rawfiles = []
# for repeat_i in range(2):
#     rawfile = f"H:\\ctzdata\\W21-8\\20240903_Ene_Calib_700fps_30mins_-700V_Am241_Co57\\20240903_Ene_Calib_700fps_30mins_-700V_Am241_Co57repeat_{repeat_i:03d}\\20240903_Ene_Calib_700fps_30mins_-700V_Am241_Co57.bin"
#     rawfiles.append(rawfile)
# output_path = r"H:\ctzdata\W21-8\processed\min30_crystals_clusters"

rawfiles = []
for repeat_i in range(3):
    rawfile = f"H:\\ctzdata\\W21-8\\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57\\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57repeat_{repeat_i:03d}\\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57.bin"
    rawfiles.append(rawfile)
output_path = r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\min45_crystals_clusters"

crystals_frames_folder = r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\crystals_frames"
crystals_spectra_folder = r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\crystals_spectra"
crystals_calibrations = [r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\calibration_crystal_0.h5"]
crystals_thresholds = [r"from_calibration"]

raw2frames(rawfiles[:1], crystals_frames_folder, config="config_czt.txt")
# raw2spectra(rawfiles, crystals_spectra_folder, config="config_czt.txt")

# raw2clusters(rawfiles, output_path, crystals_calibrations, crystals_thresholds, extended_mode=False, config="config_czt.txt")

#%% peak finding test

# n_peaks = 2
# spectra_file = r"H:\ctzdata\W21-8\processed\crystals_spectra\pixels_spectra_crystal_0.h5"
# with h5py.File(spectra_file, 'r') as f:
#     bins = f['bins'][:]
#     pixels_spectra = f['pixels_spectra'][:]

# for pixel_id in range(605, pixels_spectra.shape[0]):
#     peaks_indices = find_peaks(bins, pixels_spectra[pixel_id:pixel_id+1, :], n_peaks=2, kernel_half_length=3, neglected_bins=200, high_end_neglected_bins=200, peak_distance=20, verbose=True)
#     # peaks_indices = find_peaks(bins[:375], pixels_spectra[pixel_id:pixel_id+1, 0:375], n_peaks=1, kernel_half_length=3, neglected_bins=100, high_end_neglected_bins=0, peak_distance=5, verbose=True)
#     # peaks_indices = find_peaks(bins[375:], pixels_spectra[pixel_id:pixel_id+1, 375:], n_peaks=1, kernel_half_length=3, neglected_bins=0, high_end_neglected_bins=0, peak_distance=5, verbose=True) 
#     plt.gcf()
#     plt.title(f"Pixel {pixel_id}, Peaks at bins: {bins[peaks_indices]}")
#     # plt.pause(0.5)
#     # plt.close()
#     plt.show()

# peaks_indices = find_peaks(bins, pixels_spectra, n_peaks=2, kernel_half_length=3, neglected_bins=100, high_end_neglected_bins=200, peak_distance=20, verbose=True)
# pixels_peaks0 = peaks_indices[:, 0]
# pixels_peaks1 = peaks_indices[:, 1]


# fit_range_left = 6
# fit_range_right = 6
# for pixel_id in range(pixels_spectra.shape[0]):
#     peaks_indices = find_peaks(bins, pixels_spectra[pixel_id:pixel_id+1, :], n_peaks=n_peaks, kernel_half_length=3, neglected_bins=200, high_end_neglected_bins=200, peak_distance=5, verbose=False)
#     for peak_i in range(2):
#         peak_index = peaks_indices[0, peak_i]
#         fitted_params = fit_peaks(bins, pixels_spectra[pixel_id], peak_index, fit_range_left, fit_range_right, smooth=False, verbose=True, model='asymmetric_gaussian')
#         plt.gcf()
#         plt.title(f'Pixel {pixel_id}, Peak {peak_i}')
#         plt.legend()
#     plt.pause(0.5)
#     plt.close('all')

# fit_range_left = 6
# fit_range_right = 6
# pixels_peak0_params = np.zeros((pixels_spectra.shape[0], 7))  # amplitude, mean, sigma_left, sigma_right, background
# for pixel_id in range(pixels_spectra.shape[0]):
#     fitted_params, loss = fit_peaks(bins, pixels_spectra[pixel_id], pixels_peaks0[pixel_id], fit_range_left, fit_range_right, smooth=False, verbose=False, model='asymmetric_gaussian')
#     pixels_peak0_params[pixel_id, :6] = fitted_params
#     pixels_peak0_params[pixel_id, 6] = loss
# pixels_peak1_params = np.zeros((pixels_spectra.shape[0], 7))  # amplitude, mean, sigma_left, sigma_right, background
# for pixel_id in range(pixels_spectra.shape[0]):
#     fitted_params, loss = fit_peaks(bins, pixels_spectra[pixel_id], pixels_peaks1[pixel_id], fit_range_left, fit_range_right, smooth=False, verbose=False, model='asymmetric_gaussian')
#     pixels_peak1_params[pixel_id, :6] = fitted_params
#     pixels_peak1_params[pixel_id, 6] = loss

# with h5py.File(r"H:\ctzdata\W21-8\processed\pixels_peak0_crystal_0.h5", 'w') as f:
#     f['peak0'] = pixels_peak0_params
# with h5py.File(r"H:\ctzdata\W21-8\processed\pixels_peak1_crystal_0.h5", 'w') as f:
#     f['peak1'] = pixels_peak1_params

#%% calibration

# spectra_file = r"H:\ctzdata\W21-8\processed\crystals_spectra\pixels_spectra_crystal_0.h5"
# with h5py.File(spectra_file, 'r') as f:
#     bins = f['bins'][:]
#     pixels_spectra = f['pixels_spectra'][:]

# with h5py.File(r"H:\ctzdata\W21-8\processed\pixels_peak0_crystal_0.h5", 'r') as f:
#     pixels_peak0_params = f['peak0'][:]
# with h5py.File(r"H:\ctzdata\W21-8\processed\pixels_peak1_crystal_0.h5", 'r') as f:
#     pixels_peak1_params = f['peak1'][:]

# fig = plt.figure(figsize=(8, 6))
# ax = plt.imshow(pixels_peak0_params[:,1].reshape(40,40), cmap='gray')
# plt.colorbar(label='Peak 0 Mean (ADU)')
# plt.title('Peak 0 Mean Map')
# from pyscripts.lib.interactive_calibration_map import show_fitting_plot
# fig.canvas.mpl_connect('button_press_event', lambda event: show_fitting_plot(event, bins, pixels_spectra, asymmetric_gaussian, pixels_peak0_params, pixels_peak1_params))
# plt.show()

# pixels_start_adus = np.argmax(pixels_spectra[:, 0:], axis=1)
# # pixels_peaks1 = find_peaks(bins, pixels_spectra[:, :375], n_peaks=1, kernel_half_length=3, neglected_bins=100, high_end_neglected_bins=0, peak_distance=5, verbose=False)
# # pixels_peaks2 = find_peaks(bins, pixels_spectra[:, 375:], n_peaks=1, kernel_half_length=3, neglected_bins=0, high_end_neglected_bins=0, peak_distance=5, verbose=False) + 375
# pixels_peaks0 = pixels_peak0_params[:, 0]
# pixels_peaks1 = pixels_peak1_params[:, 0]
# energy1 = 59.54  # keV for Am241
# energy2 = 122.06  # keV for Co57

# slope, intercept, isvalids, threshold, R2 = calibrate_peaks(pixels_peaks0, pixels_peaks1, energy1, energy2, \
#         threshold_energy=5, peaks1_start_adus=pixels_start_adus, peaks2_start_adus=pixels_start_adus, verbose=True)

# calibration_filepath = r"H:\ctzdata\W21-8\processed\calibration_crystal_0.h5"
# with h5py.File(calibration_filepath, 'w') as f:
#     f.create_dataset('pixels_calibrations', data=np.vstack((slope, intercept)).T)
#     f.create_dataset('pixels_isvalids', data=isvalids)
#     f.create_dataset('pixels_peaks', data=np.vstack((pixels_peaks0, pixels_peaks1)).T)
#     f.create_dataset('pixels_thresholds', data=threshold)
#     f.create_dataset('pixels_R2', data=R2)

# from pyscripts.lib.interactive_map import InteractiveMap
# from pyscripts.lib.CrystalChecker_lib import InteractiveCrystalChecker


# interactive = InteractiveCrystalChecker("any", 1, [40, 40], bins, n_maps=5, scatter_folder=None)
# interactive.spectrum_ymax = 500
# interactive.add_subplot(0, {"Peak1": pixels_peaks0, "Peak2": pixels_peaks1, "Gain": slope.reshape(40,40), "Offset": intercept.reshape(40,40), 'Isvalid': isvalids.reshape(40,40)}, multiple_spectra={"Am241": pixels_spectra, "Co57": pixels_spectra}, \
#                          multiple_peaks={"Am241": pixels_peaks0, "Co57": pixels_peaks1}, multiple_energies={"combined": [energy1, energy2]})
# interactive.view_spectrum_in_energy = False
# interactive.if_view_scatters = False
# interactive.if_view_calibration = False
# plt.show()

# spectra_file = r"H:\ctzdata\W21-8\processed\crystals_spectra\pixels_spectra_crystal_0.h5"
# with h5py.File(spectra_file, 'r') as f:
#     bins = f['bins'][:]
#     pixels_spectra = f['pixels_spectra'][:]

# calibration_file = r"H:\ctzdata\W21-8\processed\calibration_crystal_0.h5"
# with h5py.File(calibration_file, 'r') as f:
#     slopes = f['pixels_calibrations'][:, 0]
#     intercepts = f['pixels_calibrations'][:, 1]
#     thresholds = f['pixels_thresholds'][:]
#     R2 = f['pixels_R2'][:]
#     isvalids = f['pixels_isvalids'][:]


# for row in range(40):
#     for col in range(40):
#         pixel_id = row * 40 + col
#         if not isvalids[pixel_id]:
#             continue
#         erg_bins = bins * slopes[pixel_id] + intercepts[pixel_id]
#         plt.plot(erg_bins, pixels_spectra[pixel_id])
#     plt.xlim(0, 200)
#     plt.ylim(0, 1000)
#     plt.xlabel('Energy (keV)')
#     plt.title('Calibrated Spectra for Each Pixel')
#     plt.show()

#%% thresholds
# from pyscripts.lib.calculate_thresholds import calculate_thresholds

# thresholds = calculate_thresholds(bins, pixels_spectra, percentile=0.99, verbose=False)
# with h5py.File(r"H:\ctzdata\W21-8\processed\pixels_thresholds_crystal_0.h5", 'w') as f:
#     f['pixels_thresholds'] = thresholds
