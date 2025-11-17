import h5py
import numpy as np
import matplotlib.pyplot as plt

spectra_filepath = r"H:\ctzdata\W21-8\20240903_Ene_Calib_700fps_45mins_-700V_Am241_Co57_results\crystals_spectra\pixels_spectra_crystal_0.h5"
with h5py.File(spectra_filepath, 'r') as f:
    bins = f['bins'][:]
    pixels_spectra = f['pixels_spectra'][:]

# generate the merged spectrum
merged_spectrum = np.sum(pixels_spectra, axis=0)
plt.plot(bins, merged_spectrum)
plt.yscale('log')
plt.title("Merged Spectrum")
plt.xlabel("ADU")
plt.ylabel("Counts")
plt.show()

# display a map of counts larger than 1000 ADU
counts_map = np.sum(pixels_spectra[:, bins > 1000], axis=1)
plt.imshow(counts_map.reshape(40, 40), cmap='hot')
plt.colorbar(label='Counts > 1000 ADU')
plt.title("Counts Map for ADU > 1000")
plt.show()



# display spectra pixel by pixel
for row in range(1, 40):
    for col in range(40):
        pixel_id = row * 40 + col
        plt.plot(bins, pixels_spectra[pixel_id, :])
        # plt.yscale('log')
        plt.ylim(0, 500)
        plt.title(f"Pixel ID: {pixel_id}")
        plt.xlabel("ADU")
        plt.ylabel("Counts")
        plt.show()
