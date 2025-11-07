import h5py
import numpy as np
import scipy.io as sio
import pathlib
import re
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from pyscripts.CrystalChecker_lib import InteractiveCrystalChecker

#%% Alpha spect mini with old calibration data

scatter_folder_pattern = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_{panel_id}\scatters"
calibration_file_pattern = r"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\calibration_crystal_{crystal_id}.h5"
threshold_file_pattern = r"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\pixels_thresholds_crystal_{crystal_id}.h5"
spectrum_file_Co57_pattern = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_10hrs\panel_{panel_id}\pixels_spectra_crystal_{crystal_id}.h5"
spectrum_file_Am241_pattern = r"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_4hours\panel_{panel_id}\pixels_spectra_crystal_{crystal_id}.h5"
spectrum_file_nosource_pattern = r"H:\alpha_spect_mini\nosource_500fps_3hours\panel_{panel_id}\pixels_spectra_crystal_{crystal_id}.h5"

# scatter_folder_pattern = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_{panel_id}\scatters"
# calibration_file_pattern = r"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\calibration_pedestal_crystal_{crystal_id}.h5"
# threshold_file_pattern = r"H:\alpha_spect_mini\calibration_results\panel_{panel_id}\pixels_thresholds_pedestal_crystal_{crystal_id}.h5"
# spectrum_file_Co57_pattern = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_{panel_id}\pixels_spectra_pedestal_crystal_{crystal_id}.h5"
# spectrum_file_Am241_pattern = r"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_2hours\panel_{panel_id}\pixels_spectra_pedestal_crystal_{crystal_id}.h5"
# spectrum_file_nosource_pattern = r"H:\alpha_spect_mini\nosource_500fps_3hours\panel_{panel_id}\pixels_spectra_pedestal_crystal_{crystal_id}.h5"


for panel_id in [1]:

    bins = np.linspace(6000, 14000, 1001)[:-1]
    interactive = InteractiveCrystalChecker(f"Alpha SPECT mini-Panel {panel_id}", 4, (80, 80), bins, n_maps=6)
    interactive.scatter_folder = scatter_folder_pattern.format(panel_id=panel_id)
    interactive.view_spectrum_in_energy = True
    interactive.if_view_scatters = False
    interactive.view_spectrum_in_energy = True
    interactive.if_view_calibration = False
    # interactive.spectrum_ymax = 2E6

    for crystal_id in range(4):
        # load calibration
        calibration_file = calibration_file_pattern.format(panel_id=panel_id, crystal_id=crystal_id)
        with h5py.File(calibration_file, 'r') as f:
            slopes = f['pixels_calibrations'][:, 0]
            intercepts = f['pixels_calibrations'][:, 1]
            thresholds = f['pixels_thresholds'][:]
            R2 = f['pixels_calibrations'][:, 2]
            isvalids = f['pixels_isvalids'][:]
            peaks = f['pixels_peaks'][:]

        # load thresholds
        threshold_file = threshold_file_pattern.format(panel_id=panel_id, crystal_id=crystal_id)
        with h5py.File(threshold_file, 'r') as f:
            thresholds = f['pixels_thresholds'][:]

        # load spectra
        spectrum_file_Co57 = spectrum_file_Co57_pattern.format(panel_id=panel_id, crystal_id=crystal_id)
        spectrum_file_Am241 = spectrum_file_Am241_pattern.format(panel_id=panel_id, crystal_id=crystal_id)
        spectrum_file_nosource = spectrum_file_nosource_pattern.format(panel_id=panel_id, crystal_id=crystal_id)
        with h5py.File(spectrum_file_Co57, 'r') as f:
            spectra_Co57 = f['pixels_spectra'][:]
        with h5py.File(spectrum_file_Am241, 'r') as f:
            spectra_Am241 = f['pixels_spectra'][:]
        with h5py.File(spectrum_file_nosource, 'r') as f:
            spectra_nosource = f['pixels_spectra'][:]
        multiple_spectra = {"Co57": spectra_Co57, "Am241": spectra_Am241, 'nosource': spectra_nosource}
        multiple_peaks = {"Co57": peaks[:,  1:], "Am241": peaks[:, 0], "nosource": np.zeros((6400, 1), dtype=np.float32)}
        multiple_energies = {"Co57": [122.0615], "Am241": 59.5409, "nosource": 0.0}

        # load hasscatters
        hasscatters = np.zeros((80, 80), dtype=bool)
        files = [file.name for file in pathlib.Path(interactive.scatter_folder).glob("scatter_crystal_*.h5")]
        crystal_ids = [int(re.match(r"scatter_crystal_(\d+)_pixel_(\d+).h5", file).group(1)) for file in files]
        pixel_ids =  [int(re.match(r"scatter_crystal_(\d+)_pixel_(\d+).h5", file).group(2)) for file in files]
        for icrystal_id, pixel_id in zip(crystal_ids, pixel_ids):
            if icrystal_id == crystal_id:
                row = pixel_id // 80
                col = pixel_id % 80
                hasscatters[row, col] = 1

        # add subplot
        interactive.add_subplot(crystal_id, {"Gain": slopes, "Offset": intercepts, "R2": R2, "Thresholds": thresholds, "IsValid": isvalids, "HasScatters": hasscatters},
                                 multiple_peaks, multiple_spectra, multiple_energies)
    plt.show()



