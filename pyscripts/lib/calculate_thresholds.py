import numpy as np
import matplotlib.pyplot as plt


def calculate_thresholds(bins, pixels_spectra, percentile=0.99, verbose=False):
    """
    Calculate the threshold for each pixel based on the given percentile of its spectrum.

    Parameters:
    spectra (np.ndarray): 2D array where each row corresponds to a pixel's spectrum.
    bins (np.ndarray): 1D array of bin edges corresponding to the spectra.
    percentile (float): The percentile to use for threshold calculation.

    Returns:
    np.ndarray: 1D array of thresholds for each pixel.
    """
    pixels_cumsum = np.cumsum(pixels_spectra, axis=1)
    pixels_total_counts = pixels_cumsum[:, -1]  # 每个像素的总和
    targets = pixels_total_counts * (percentile)
    pixels_threshold_indices = np.argmax(pixels_cumsum >= targets[:, None], axis=1)
    pixels_higher_counts = pixels_cumsum[np.arange(pixels_cumsum.shape[0]), pixels_threshold_indices-1]
    pixels_remaining_counts = targets - pixels_higher_counts
    pixels_thresholds = bins[pixels_threshold_indices] + np.round(pixels_remaining_counts / pixels_spectra[np.arange(pixels_spectra.shape[0]), pixels_threshold_indices] * 8)

    if verbose:
        for pixel_i in range(min(10, pixels_spectra.shape[0])):
            print(f"Pixel {pixel_i} - Total Counts: {pixels_total_counts[pixel_i]}({percentile*100:.1f}% ={percentile*pixels_total_counts[pixel_i]}):\
{pixels_thresholds[pixel_i]}, {sum(pixels_spectra[pixel_i, :pixels_threshold_indices[pixel_i]])} counts lower than threshold,\
{sum(pixels_spectra[pixel_i, :pixels_threshold_indices[pixel_i]]) / pixels_total_counts[pixel_i] * 100:.1f}%")

            plt.figure(figsize=(8, 5))
            plt.plot(bins, pixels_spectra[pixel_i, :]+1, label=f'Pixel {pixel_i} Spectrum', marker='o')
            plt.axvline(pixels_thresholds[pixel_i], color='red', linestyle='--', label=f'{percentile}th Percentile Threshold')
            plt.title(f'Pixel {pixel_i} Spectrum with Threshold')
            plt.xlabel('Bins')
            plt.ylabel('Counts')
            plt.xlim(pixels_thresholds[pixel_i]-200, pixels_thresholds[pixel_i]+200)
            plt.legend()
            plt.show()

    return pixels_thresholds


if __name__ == "__main__":

    import h5py

    # Example usage
    spectrum_file = r"H:\alpha_spect_mini\nosource_500fps_3hours\panel_1\pixels_spectra_crystal_0.h5"
    with h5py.File(spectrum_file, 'r') as f:
        bins = f['bins'][:]
        pixels_spectra = f['pixels_spectra'][:]

    pixels_thresholds = calculate_thresholds(bins, pixels_spectra, percentile=0.99, verbose=True)
