import pathlib
from pyscripts.lib.cmd_wrapper import *

panel_id = 1

# rawfolder = pathlib.Path(rf"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_{panel_id}")
# rawfiles = list(rawfolder.glob(f"20241120_EneCalib_Co57_150fps_-500V_7hrs_*_0{panel_id-1}.bin"))
# rawfiles = sorted([str(f) for f in rawfiles])


# rawfolder = pathlib.Path(rf"H:\alpha_spect_mini\20250218_pc123456_usb123456_daq123456_30min_Am241\panel_{panel_id}")
# rawfiles = list(rawfolder.glob(f"test_0218_pc123456_usb123456_daq123456_*_0{panel_id-1}.bin"))
# rawfiles = sorted([str(f) for f in rawfiles])

# rawfolder = pathlib.Path(rf"H:\alpha_spect_mini\20250225_pc123456_usb123456_daq123456_800fps_30min_Co57\panel_{panel_id}")
# rawfiles = list(rawfolder.glob(f"20250225_Co57_60min_800fps_*_0{panel_id-1}.bin"))
# rawfiles = sorted([str(f) for f in rawfiles])

rawfolder = pathlib.Path(rf"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_2hours\panel_{panel_id}")
rawfiles = list(rawfolder.glob(f"20241220_Am241_panel12356_800fps_2hours_*_0{panel_id-1}.bin"))
rawfiles = sorted([str(f) for f in rawfiles])


crystals_calibrations = pathlib.Path(rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}").glob("calibration_crystal_?.h5")
crystals_calibrations = sorted([str(f) for f in crystals_calibrations])
# crystals_calibrations[0] = 'skip'

crystals_thresholds = pathlib.Path(rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}").glob("pixels_thresholds_crystal_?.h5")
crystals_thresholds = sorted([str(f) for f in crystals_thresholds])
# crystals_thresholds = ['from_calibration' for _ in crystals_thresholds]

crystals_spectra_folder = rawfolder / "crystals_spectra"
crystals_spectra_folder = str(crystals_spectra_folder)

crystals_cluster_folder = rawfolder / "crystals_clusters"
crystals_cluster_folder = str(crystals_cluster_folder)

crystals_frames_folder = rawfolder / "crystals_frames"
crystals_frames_folder = str(crystals_frames_folder)

crystals_scatters_folder = rawfolder / "crystals_scatters"
crystals_scatters_folder = str(crystals_scatters_folder)

# scattering_mode = "all-stride"
extended_mode = False

# raw2frames(rawfiles[:1], crystals_frames_folder, cwd=r"F:\alpha\cpp_plugins")
# raw2spectra(rawfiles[:1], crystals_spectra_folder, cwd=r"F:\alpha\matlab\select_clusters\cpp_plugins")
# raw2scatters(rawfiles[:5], crystals_scatters_folder, crystals_calibrations, mode=scattering_mode, config="config_alpha")
raw2clusters(rawfiles, crystals_cluster_folder, crystals_calibrations, crystals_thresholds, extended_mode, config="config_alpha.txt")