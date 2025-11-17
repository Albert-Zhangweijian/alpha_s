from pyscripts.lib.read_clusters import read_clusters, write_clusters

Co57_calibration_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\Co57_pixels_calibration_crystal_0.h5"
import h5py
with h5py.File(Co57_calibration_filepath, 'r') as f:
    pixels_calibration = f['pixels_calibrations'][:]


clusters = read_clusters(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\clusters_crystal_0.bin")


for cluster in clusters:
    for i in range(len(cluster.energies)):
        pixel_id = cluster.pixel_ids[i]
        adu = cluster.adus[i]
        energy = adu * pixels_calibration[pixel_id, 0] + pixels_calibration[pixel_id, 1]
        cluster.energies[i] = energy  # convert to keV

write_clusters(r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\clusters_crystal_0_Co57_reset_energies.bin", clusters)