import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve1d

def generate_centroid_kernel(kernel_half_length):
    if kernel_half_length == 0:
        return np.array([1], dtype=float)
    else:
        previous_kernel = generate_centroid_kernel(kernel_half_length - 1)
        kernel = np.zeros(2 * kernel_half_length + 1, dtype=float)
        kernel[:-2] = previous_kernel / 4
        kernel[1:-1] += previous_kernel / 2
        kernel[2:] += previous_kernel / 4
        return kernel

def generate_0d_savitzky_kernel(kernel_half_length):
    if kernel_half_length == 1:
        return np.array([0, 1, 0], dtype=float)
    if kernel_half_length == 2:
        return np.array([-0.085714, 0.342857, 0.485714, 0.342857, -0.085714], dtype=float)
    if kernel_half_length == 3:
        return np.array([-0.095238, 0.142857, 0.285714, 0.333333, 0.285714, 0.142857, -0.095238], dtype=float)
    if kernel_half_length == 4:
        return np.array([-0.090909, 0.060606, 0.168831, 0.233766, 0.255411, 0.233766, 0.168831, 0.060606, -0.090909], dtype=float)
    if kernel_half_length == 5:
        return np.array([-0.083916, 0.020979, 0.102564, 0.160839, 0.195804, 0.207459, 0.195804, 0.160839, 0.102564, 0.020979, -0.083916], dtype=float)
    if kernel_half_length == 6:
        return np.array([-0.076923, 0.000000, 0.062937, 0.111888, 0.146853, 0.167832, 0.174825, 0.167832, 0.146853, 0.111888, 0.062937, 0.000000, -0.076923], dtype=float)
    if kernel_half_length == 7:
        return np.array([-0.070588, -0.011765, 0.038009, 0.078733, 0.110407, 0.133032, 0.146606, 0.151131, 0.146606, 0.133032, 0.110407, 0.078733, 0.038009, -0.011765, -0.070588], dtype=float)
    raise ValueError("Unsupported kernel_half_length(not tabulated): {}".format(kernel_half_length))


def generate_1d_savitzky_kernel(kernel_half_length):
    if kernel_half_length == 1:
        return np.array([-0.5, 0, 0.5], dtype=float)
    if kernel_half_length == 2:
        return np.array([-0.200000, -0.100000, 0.000000, 0.100000, 0.200000], dtype=float)
    if kernel_half_length == 3:
        return np.array([-0.107143, -0.071429, -0.035714, 0.0, 0.035714, 0.071429, 0.107143], dtype=float)
    if kernel_half_length == 4:
        return np.array([-0.06666667, -0.05, -0.03333333, -0.01666667, 0.0, 0.01666667, 0.03333333, 0.05, 0.06666667], dtype=float)
    if kernel_half_length == 5:
        return np.array([-0.045455, -0.036364, -0.027273, -0.018182, -0.009091, 0.0, 0.009091, 0.018182, 0.027273, 0.036364, 0.045455], dtype=float)
    if kernel_half_length == 6:
        return np.array([-0.032967, -0.027473, -0.021978, -0.016484, -0.010989, -0.005495, 0.0, 0.005495, 0.010989, 0.016484, 0.021978, 0.027473, 0.032967], dtype=float)
    if kernel_half_length == 7:
        return np.array([-0.025000, -0.021429, -0.017857, -0.014286, -0.010714, -0.007143, -0.003571, 0.0, 0.003571, 0.007143, 0.010714, 0.014286, 0.017857, 0.021429, 0.025000], dtype=float)
    raise ValueError("Unsupported kernel_half_length(not tabulated): {}".format(kernel_half_length))

def find_peaks(bins, spectrum, n_peaks, smooth_kernel_half_length=3, savitzky_kernel_half_length=3, neglected_bins=200, high_end_neglected_bins=200, peak_distance=5, verbose=False):

    assert n_peaks > 0
    assert peak_distance > 0

    peaks_indices = np.zeros((spectrum.shape[0], n_peaks), dtype=int)
    centroid_kernel = generate_centroid_kernel(smooth_kernel_half_length)
    savitzky_kernel = generate_1d_savitzky_kernel(savitzky_kernel_half_length)
    ymax = np.max(spectrum[:, neglected_bins:])  # to ignore the low-energy noise peak

    # smooth the spectrum with Savitzky-Golay filter
    smoothed_spectrum = convolve1d(spectrum, centroid_kernel, axis=-1, mode='constant')
    convolved_spectrum = convolve1d(smoothed_spectrum, savitzky_kernel, axis=-1, mode='constant')

    # check_start_index = np.zeros(spectrum.shape[0], dtype=int) + neglected_bins
    # for peak_i in range(n_peaks):
    #     convolved_spectrum_copy = np.copy(convolved_spectrum)
    #     for pixel_i in range(spectrum.shape[0]):
    #         convolved_spectrum_copy[pixel_i, 0:check_start_index[pixel_i]] = 1E5  # mask the area before check_start_index
    #     peak_index = np.argmin(convolved_spectrum_copy, axis=1)
    #     peaks_indices[:, peak_i] = peak_index + kernel_half_length  # adjust for the kernel shift
    #     check_start_index = peak_index + peak_distance
    #     check_start_index[peak_index+peak_distance >= (spectrum.shape[1] - high_end_neglected_bins)] = len(spectrum) - high_end_neglected_bins - 1

    convolved_spectrum_copy = np.copy(convolved_spectrum)
    convolved_spectrum_copy[:, :neglected_bins] = 0  # mask the area before neglected_bins
    convolved_spectrum_copy[:, spectrum.shape[1]-high_end_neglected_bins:] = 0  # mask the area after high_end_neglected_bins

    for peak_i in range(n_peaks):
        peak_index = np.argmax(convolved_spectrum_copy, axis=1)
        peaks_indices[:, peak_i] = peak_index - savitzky_kernel_half_length / 2 # adjust for the kernel shift
        left_mask_indices = peak_index - peak_distance
        left_mask_indices[left_mask_indices < neglected_bins] = neglected_bins
        right_mask_indices = peak_index + peak_distance + 1
        right_mask_indices[right_mask_indices > spectrum.shape[1] - high_end_neglected_bins] = spectrum.shape[1] - high_end_neglected_bins
        for i in range(spectrum.shape[0]):
            convolved_spectrum_copy[i, left_mask_indices[i]:right_mask_indices[i]] = 0
    peaks_indices = np.sort(peaks_indices, axis=1)


    if verbose and spectrum.shape[0] == 1:
        plt.figure(figsize=(10, 6))
        plt.plot(bins,np.squeeze(spectrum), label='Spectrum', linestyle='none', marker='.', color='blue', markersize=4)
        plt.plot(bins, np.squeeze(smoothed_spectrum), label='Smoothed Spectrum', linestyle='-', color='orange')
        plt.plot(bins, np.squeeze(convolved_spectrum), label='Convolved Spectrum', linestyle='-', color='green')
        plt.scatter(bins[peaks_indices], spectrum[0, peaks_indices], color='red', label='Detected Peaks', marker='d')
        plt.axvline(x=bins[neglected_bins], color='green', linestyle='--', label='Neglected Bins Start')
        plt.axvline(x=bins[len(spectrum) - high_end_neglected_bins], color='green', linestyle='--', label='Neglected Bins End')
        plt.title(f'Spectrum with Detected Peaks')
        plt.xlabel('Bins')
        plt.ylabel('Counts')
        plt.ylim(0, ymax*1.1)
        plt.legend()
        plt.grid()
        plt.tight_layout()

    return peaks_indices


if __name__ == "__main__":

    import h5py
    

    n_peaks = 1
    spectra_file = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_spectra\pixels_spectra_crystal_0.h5"
    with h5py.File(spectra_file, 'r') as f:
        bins = f['bins'][:]
        pixels_spectra = f['pixels_spectra'][:]

    for pixel_id in range(700, pixels_spectra.shape[0]):
        peaks_indices = find_peaks(bins, pixels_spectra[pixel_id:pixel_id+1, :], n_peaks=n_peaks, kernel_half_length=3, neglected_bins=200, high_end_neglected_bins=200, peak_distance=5, verbose=True)
        plt.gcf()
        plt.title(f"Pixel {pixel_id}, Peaks at bins: {bins[peaks_indices]}")
        # plt.pause(1)
        # plt.close()
        plt.show()
