import h5py
import numpy as np
import matplotlib.pyplot as plt

from lib.read_clusters import read_clusters_fast

# frames_filepath = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_frames\frames_0_40000_crystal_0.h5"
# calibrations_filepath = r"H:\alpha_spect_mini\calibration_results\panel_1\calibration_crystal_0_new.h5"
# clusters_folderpath = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_clusters"

frames_filepath = r"H:\ctzdata\W21-8\processed\crystals_frames\frames_0_10000_crystal_0.h5"
calibrations_filepath = r"H:\ctzdata\W21-8\processed\calibration_crystal_0.h5"
clusters_folderpath = r"H:\ctzdata\W21-8\processed\min45_crystals_clusters"
shape = [40, 40]

with h5py.File(frames_filepath, 'r') as f:
    frames = f['frames'][:]

with h5py.File(calibrations_filepath, 'r') as f:
    slopes = f['pixels_calibrations'][:, 0]
    intercepts = f['pixels_calibrations'][:, 1]
    isvalids = f['pixels_isvalids'][:]
erg_frames = frames * slopes.T + intercepts.T  # Convert to energy using calibration of crystal 1

clusters_npixels_dict = {}
for n_pixels in range(1, 10):
    clusters_filepath = clusters_folderpath + rf"\clusters_crystal_0_pixel_{n_pixels}.bin"
    clusters_npixels_dict[n_pixels] = read_clusters_fast(clusters_filepath)

bad_pixels = np.where(~isvalids)[0]
bad_pixels_rows = bad_pixels // shape[1]
bad_pixels_cols = bad_pixels % shape[1]

n_frames = frames.shape[0]
for frame_id in range(0, n_frames):

    frame = frames[frame_id].reshape(shape)
    erg_frame = erg_frames[frame_id].reshape(shape)
    plt.figure(figsize=(8, 8))
    plt.imshow(erg_frame, cmap='gray', vmin=0, vmax=200)

    plt.scatter(bad_pixels_cols, bad_pixels_rows, marker='x', color='blue', label='Bad Pixels')

    # Get clusters corresponding to this frame
    for n_pixels in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        clusters_dict = clusters_npixels_dict[n_pixels]
        first_file_ending_index = np.where(clusters_dict['frame_ids'][1:] < clusters_dict['frame_ids'][:-1])[0][0]

        cluster_indices = np.where(clusters_dict['frame_ids'] == frame_id)[0]
        cluster_indices = cluster_indices[cluster_indices < first_file_ending_index]
        pixel_ids = clusters_dict['pixel_ids'][:, cluster_indices]
        adus = clusters_dict['adus'][:, cluster_indices]

        for i in range(len(cluster_indices)):
            x_coords = pixel_ids[:, i] % shape[1]
            y_coords = pixel_ids[:, i] // shape[1]
            x_edges = [x_coords - 0.5, x_coords + 0.5, x_coords + 0.5, x_coords - 0.5, x_coords - 0.5]
            y_edges = [y_coords - 0.5, y_coords - 0.5, y_coords + 0.5, y_coords + 0.5, y_coords - 0.5]
            plt.plot(x_edges, y_edges, color='red')

    plt.title(f"Frame ID: {frame_id}")
    plt.colorbar(label='ADU')
    plt.show()