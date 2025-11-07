# import numpy as np
# import h5py
# from find_peaks import find_peaks
# from fit_peaks import fit_peaks


# for panel_i in range(6):
#     for crystal_i in range(4):
#         spectra = f"H:\\alpha_spect_mini\\20241120_EneCalib_Co57_150fps_-500V_7hrs\\panel_{panel_i}\\pixels_spectra_crystal_{crystal_i}.h5"

#         with h5py.File(spectra, 'r') as f:
#             bins = f['bins'][:]
#             pixels_spectra = f['pixels_spectra'][:]
#             n_pixels = pixels_spectra.shape[0]

#             for pixel_id in range(n_pixels):
#                 peaks_indices = find_peaks(bins, pixels_spectra[pixel_id:pixel_id+1, :], n_peaks=n_peaks, kernel_half_length=3, neglected_bins=200, high_end_neglected_bins=200, peak_distance=5, verbose=True)
#                 plt.gcf()
#                 plt.title(f"Pixel {pixel_id}, Peaks at bins: {bins[peaks_indices]}")
#                 # plt.pause(1)
#                 # plt.close()
#                 plt.show()
