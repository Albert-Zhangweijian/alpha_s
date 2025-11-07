import pathlib
from pyscripts.cmd_wrapper import *

panel_id = 1

rawfolder = pathlib.Path(rf"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_{panel_id}")
rawfiles = list(rawfolder.glob(f"20241120_EneCalib_Co57_150fps_-500V_7hrs_*_0{panel_id-1}.bin"))
rawfiles = sorted([str(f) for f in rawfiles])

crystals_calibrations = pathlib.Path(rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}").glob("calibration_crystal_?.h5")
crystals_calibrations = sorted([str(f) for f in crystals_calibrations])

crystals_thresholds = pathlib.Path(rf"H:\alpha_spect_mini\calibration_results\panel_{panel_id}").glob("pixels_thresholds_crystal_?.h5")
crystals_thresholds = sorted([str(f) for f in crystals_thresholds])

crystals_spectra_folder = rawfolder / "crystals_spectra"
crystals_spectra_folder = str(crystals_spectra_folder)

crystals_cluster_folder = rawfolder / "crystals_clusters"
crystals_cluster_folder = str(crystals_cluster_folder)

crystals_frames_folder = rawfolder / "crystals_frames"
crystals_frames_folder = str(crystals_frames_folder)

crystals_scatters_folder = rawfolder / "crystals_scatters"
crystals_scatters_folder = str(crystals_scatters_folder)


clustering_mode = "random"

# raw2frames(rawfiles[:1], crystals_frames_folder, cwd=r"F:\alpha\matlab\select_clusters\cpp_plugins")
# raw2spectra(rawfiles[:1], crystals_spectra_folder, cwd=r"F:\alpha\matlab\select_clusters\cpp_plugins")
# raw2scatters(rawfiles[:1], crystals_scatters_folder, crystals_calibrations, mode='all-random', cwd=r"F:\alpha\matlab\select_clusters\cpp_plugins")
raw2clusters(rawfiles, crystals_cluster_folder, crystals_calibrations, crystals_thresholds, cwd=r"F:\alpha\matlab\select_clusters\cpp_plugins")