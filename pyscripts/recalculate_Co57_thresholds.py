import numpy as np
import matplotlib.pyplot as plt
import h5py


filepath = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_spectra\pixels_spectra_crystal_0.h5"
with h5py.File(filepath, 'r') as f:
    bins = f['bins'][:]
    pixels_spectra = f['pixels_spectra'][:]
    
calibration_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\calibration_crystal_0.h5"
with h5py.File(calibration_filepath, 'r') as f:
    pixels_calibration = f['pixels_calibrations'][:]
    pixels_peaks = f['pixels_peaks'][:]


first_nonzero_indices = np.argmax(pixels_spectra > 0, axis=1)
first_nonzero_adus = bins[first_nonzero_indices] + 4

photopeak_adus = pixels_peaks[:, 1]

pixels_gain = 122.06 / (photopeak_adus - first_nonzero_adus)  # eV/ADU
pixels_offset = -pixels_gain * first_nonzero_adus  # eV
pixels_thresholds = (5 - pixels_offset) / pixels_gain  # ADU

counts_lower_than_5kev = np.array([np.sum(spectrum[bins * pixels_gain[pixel_id] + pixels_offset[pixel_id] <= 5]) for pixel_id, spectrum in enumerate(pixels_spectra)])
counts_lower_than_10kev = np.array([np.sum(spectrum[bins * pixels_gain[pixel_id] + pixels_offset[pixel_id] <= 10]) for pixel_id, spectrum in enumerate(pixels_spectra)])

new_pixels_calibration = np.vstack((pixels_gain, pixels_offset)).T

# Co57_calibration_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\Co57_pixels_calibration_crystal_0.h5"
# with h5py.File(Co57_calibration_filepath, 'w') as f:
#     f.create_dataset('pixels_calibrations', data=new_pixels_calibration)



for pixel_id in range(pixels_spectra.shape[0]):
    erg_bins = bins * new_pixels_calibration[pixel_id, 0] + new_pixels_calibration[pixel_id, 1]
    print(f"Counts lower than 5 eV per pixel is {counts_lower_than_5kev[pixel_id]}, ratio is {counts_lower_than_5kev[pixel_id] / counts_lower_than_10kev[pixel_id]}")
    
    plt.figure(figsize=(10, 6))
    plt.plot(erg_bins, pixels_spectra[pixel_id], label=f'Pixel {pixel_id}')
    plt.axvline(5, color='r', linestyle='--', label='5 eV Threshold')
    plt.xlim(0, 200)  # Adjust x-axis limit as needed
    plt.ylim(0, 1000)
    plt.xlabel('Energy (eV)')
    plt.ylabel('Counts')
    plt.title('Pixel Spectra with New Calibration')
    plt.legend()
    plt.show()
