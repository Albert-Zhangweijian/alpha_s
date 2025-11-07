import numpy as np
import matplotlib.pyplot as plt

import h5py
from interactive_map import *

"""
This script is to check the percentage of counts below the threshold for each pixel.
"""

if_plot_single_pixel_spectrum = 1
percentile = 0.99

for panel in range(1, 7):
    for crystal in range(4):

        thresholds_file = rf"H:\alpha_spect_mini\calibration_results\panel_{panel}\pixels_thresholds_crystal_{crystal}.h5"
        spectra_file = rf"H:\alpha_spect_mini\nosource_500fps_3hours\panel_{panel}\pixels_spectra_crystal_{crystal}.h5"
    
        with h5py.File(thresholds_file, "r") as f:
            thresholds = np.array(f['pixels_thresholds'])
        with h5py.File(spectra_file, "r") as f:
            pixels_spectra = np.array(f['pixels_spectra'])
            bins = np.array(f['bins'])        

        pixels_thresholds_ratios = []
        for pixel_id in range(6400):

            # calculate the ratio of counts below the threshold
            accumulated_counts = 0
            for i in range(1000-1):
                bin_lb = bins[i]
                bin_ub = bins[i+1]
                count = pixels_spectra[pixel_id][i]
                if bin_ub <= thresholds[pixel_id]:
                    accumulated_counts += count
                elif bin_ub > thresholds[pixel_id] and bin_lb < thresholds[pixel_id]:
                    accumulated_counts += (thresholds[pixel_id]-bin_lb) / (bin_ub-bin_lb) * count
                    break
                else:
                    break

            # plot the spectrum of one pixel, and 0.99 cumulative ratio line, and calculated threshold line
            if if_plot_single_pixel_spectrum:
                total_count = pixels_spectra[pixel_id].sum()
                threshold_count = total_count * percentile
                ratios = np.cumsum(pixels_spectra[pixel_id]) / total_count
                plt.plot([bin_lb, bin_ub], [np.sum(pixels_spectra[pixel_id][:i])/total_count, np.sum(pixels_spectra[pixel_id][:i+1])/total_count], color='red', linestyle='--', label='Accumulated Counts')
                plt.step(bins, ratios, where='post', color='blue', label='Cumulative Ratio')
                plt.hlines(0.99, 6000, 14000, color='orange', linestyle='--')
                plt.vlines(thresholds[pixel_id], 0, 1, color='orange', linestyle='--')
                plt.xlim(6000, 8000)
                plt.show()
            


            pixels_thresholds_ratios.append(accumulated_counts / sum(pixels_spectra[pixel_id]))
            if pixel_id < 40:
                print(f"Pixel id = {pixel_id:3d}, Threshold={thresholds[pixel_id]}, Counts Ratio = f{accumulated_counts / sum(pixels_spectra[pixel_id]) * 100:.2f}%")

        # first plot the histogram of the ratios in one crystal
        # plt.hist(pixels_thresholds_ratios, bins=np.linspace(0.9, 1, 100))
        # plt.show()
        
        # the cumulated ratio map
        r = ThresholdMap(np.array(pixels_thresholds_ratios), pixels_spectra, maprange=(0.98, 1), thresholds=thresholds)
        r.spectrum_range = (6000, 8000)

        plt.show()
        
