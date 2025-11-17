import h5py
import numpy as np
import matplotlib.pyplot as plt

from lib.read_clusters import read_clusters_fast

frames_filepath = r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\crystals_frames\frames_160000_170000_crystal_0.h5"

with h5py.File(frames_filepath, 'r') as f:
    frames = np.asarray(f['frames'][:])

# generate the merged spectrum
merged_spectrum = np.histogram(frames.flatten(), bins=8000, range=(0, 8000))[0]
plt.plot(merged_spectrum)
plt.yscale('log')
plt.title("Merged Spectrum")
plt.xlabel("ADU")
plt.ylabel("Counts")
plt.show()

# # display frames frame by frame
# n_frames = frames.shape[0]
# for frame_id in range(0, n_frames):

#     frame = frames[frame_id].reshape(40, 40)
#     plt.figure(figsize=(8, 8))
#     plt.imshow(frame, cmap='gray', vmin=0, vmax=200)
#     plt.title(f"Frame ID: {frame_id}")
#     plt.colorbar(label='ADU')
#     plt.show()

# display response pixel by pixel
n_pixels = frames.shape[1]
for pixel_id in range(0, n_pixels):

    pixel_scatter = frames[:, pixel_id]
    plt.figure(figsize=(8, 6))
    plt.plot(pixel_scatter, '.')
    plt.title(f"Pixel ID: {pixel_id} Scatter")
    plt.xlabel("Frame Index")
    plt.ylabel("ADU")
    plt.show()