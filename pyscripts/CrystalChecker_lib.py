import h5py
import numpy as np
import scipy.io as sio
import pathlib
import re
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#%% Definition of the InteractiveCrystalChecker class

class InteractiveCrystalChecker:

    def __init__(self, title, n_crystals, shape, bins, n_maps=5, scatter_folder=None):
        
        # macro info
        self.title = title
        self.n_crystals = n_crystals
        self.n_row, self.n_col = shape
        self.bins = bins
        self.n_maps = n_maps
        self.scatter_folder = scatter_folder

        # set the figure
        self.fig, self.axes = plt.subplots(n_crystals, self.n_maps, figsize=(4 * self.n_maps, 4 * n_crystals))
        self.fig.suptitle(title)
        self.axes = np.array(self.axes).reshape(n_crystals, self.n_maps)

        # global data containers
        self.crystals_maps = [_ for _ in range(n_crystals)]
        self.all_multiple_spectra = [_ for _ in range(n_crystals)]
        self.all_multiple_peaks = [_ for _ in range(n_crystals)]
        self.all_multiple_energies = [_ for _ in range(n_crystals)]

        # configuration
        self.if_view_spectrum = True
        self.if_view_calibration = True
        self.if_view_scatters = True
        self.if_show_threshold_in_spectrum = True
        self.if_show_peaks_in_spectrum = True
        self.if_show_threshold_in_scatters = True
        self.if_show_peaks_in_scatters = True

        self.view_spectrum_in_energy = False

        self.spectrum_ymax = 100
        self.spectrum_xrange = (-5, 200)

        # automatic params
        self.map_ranges = {"Gain": (0, 0.05), "Offset": (-400, 10), "Thresholds": (6000, 8000), "R2": (0.95, 1), "IsValids": (0, 1), "HasScatters": (0, 1)}
        self.map_cmaps = {"Gain": "hot", "Offset": "hot", "Thresholds": "hot", "R2": "hot", "IsValids": "gray", "HasScatters": "gray"}
        

    def add_subplot(self, crystal_id, maps, multiple_peaks: dict, multiple_spectra: dict, multiple_energies: dict):

        # data in dict[map_title, map]
        self.crystals_maps[crystal_id] = maps

        # data in dict[nuc_name, map]
        self.all_multiple_spectra[crystal_id] = multiple_spectra
        self.all_multiple_peaks[crystal_id] = multiple_peaks
        self.all_multiple_energies[crystal_id] = multiple_energies

        # draw maps
        titles = [title for title in maps.keys()]
        for image_i, title in enumerate(titles):

            # draw the image
            map = maps[title]
            
            cmap = self.map_cmaps[title] if title in self.map_cmaps.keys() else "hot"
            vmin = self.map_ranges[title][0] if title in self.map_ranges.keys() else None
            vmax = self.map_ranges[title][1] if title in self.map_ranges.keys() else None
            
            im = self.axes[crystal_id][image_i].imshow(map.reshape((self.n_row, self.n_col)), cmap=cmap, origin="upper", interpolation='none', vmin=vmin, vmax=vmax)
            self.fig.colorbar(im, ax=self.axes[crystal_id][image_i], fraction=0.046, pad=0.04)
            self.axes[crystal_id][image_i].set_aspect('equal')
            self.axes[crystal_id][image_i].set_title(f"Crystal {crystal_id+1}/{self.n_crystals}-{title}")

            # add gridlines
            grid_lines = []
            for i in range(self.n_row):
                grid_lines.append([(-0.5, i+0.5), (self.n_col-0.5, i+0.5)])
            for j in range(self.n_col):
                grid_lines.append([(j+0.5, -0.5), (j+0.5, self.n_row-0.5)])
            grid = LineCollection(grid_lines, color="black", linewidths=0.2)
            self.axes[crystal_id][image_i].add_collection(grid)

            # set labels
            if crystal_id == self.n_crystals - 1:
                self.axes[crystal_id][image_i].set_xlabel("Pixel X")
            if image_i == 0:
                self.axes[crystal_id][image_i].set_ylabel("Pixel Y")

        # add click event
        # for spectra, we can not pass it as a parameter to on_click, so we need to store it in the class as all_multiple_spectra, etc.
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.update_plot()

    def update_plot(self):
        self.fig.canvas.draw()

    def on_click(self, event):

        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode != '':  # when using zoom, pan, etc. ignore clicks
            print(f"Click ignored: {toolbar.mode} mode is active.")
            return

        if event.xdata is None or event.ydata is None:
            return

        row, col = round(event.ydata), round(event.xdata)
        print(f"Clicked on ({event.ydata:.2f}, {event.xdata:.2f}) => Row={row} Col={col}")
        if event.inaxes:  # check if the click is inside the axes
            for i, ax in enumerate(self.axes.flatten()):
                if event.inaxes == ax:
                    print(f"Clicked on Axes {i}, i.e. Crystal {i//4+1}, {ax.get_title()}")
                    crystal_i = i // self.n_maps
                    map_i = i % self.n_maps
                    break
        for ax in self.axes[crystal_i]:
            point = ax.plot(col, row, 'cs', fillstyle='none', markersize=4)[0]

        if self.if_view_spectrum:
            self.view_pixels_spectrum(row, col, crystal_i)

        if self.if_view_calibration:
            self.view_pixels_calibration(row, col, crystal_i)

        if self.if_view_scatters:
            self.view_pixels_scatters(row, col, crystal_i)

        self.update_plot()

    def view_pixels_spectrum(self, row, col, crystal_id):
        fig, axes = plt.subplots(3, 3, figsize=(9, 9))
        plt.suptitle(f"Crystal {crystal_id} - Spectrum")

        # draw the spectrum
        for name, spectra in self.all_multiple_spectra[crystal_id].items():
            for i_row in [-1, 0, 1]:
                for i_col in [-1, 0, 1]:
                    if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                        index_spectrum = (row + i_row) * self.n_col + (col + i_col)

                        if self.view_spectrum_in_energy:
                            energies = self.bins * self.crystals_maps[crystal_id]['Gain'][index_spectrum] + self.crystals_maps[crystal_id]['Offset'][index_spectrum]
                            axes[i_row+1][i_col+1].plot(energies, spectra[index_spectrum, :], label=name)
                            axes[i_row+1][i_col+1].set_title(f"Pixel ({row + i_row}, {col + i_col})")
                            if self.if_show_peaks_in_spectrum and name != 'nosource':
                                axes[i_row+1][i_col+1].vlines(self.all_multiple_peaks[crystal_id][name][index_spectrum] * self.crystals_maps[crystal_id]['Gain'][index_spectrum] + self.crystals_maps[crystal_id]['Offset'][index_spectrum], 0, 1200, color="red", label="Peak-"+name, linestyle='--')
                            axes[i_row+1][i_col+1].set_xlim(self.spectrum_xrange)
                        else:
                            axes[i_row+1][i_col+1].plot(self.bins, spectra[index_spectrum, :], label=name)
                            axes[i_row+1][i_col+1].set_title(f"Pixel ({row + i_row}, {col + i_col})")
                            if self.if_show_peaks_in_spectrum and name != 'nosource':
                                axes[i_row+1][i_col+1].vlines(self.all_multiple_peaks[crystal_id][name][index_spectrum], 0, 1200, linestyle='--', color="red", label="Peak-"+name)

        # plot thresholds
        for i_row in [-1, 0, 1]:
            for i_col in [-1, 0, 1]:
                if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                    index_spectrum = (row + i_row) * self.n_col + (col + i_col)
                    if self.view_spectrum_in_energy:
                        axes[i_row+1][i_col+1].set_xlabel("energy (keV)")
                        if self.if_show_threshold_in_spectrum and 'Thresholds' in self.crystals_maps[crystal_id].keys():
                            axes[i_row+1][i_col+1].vlines(self.crystals_maps[crystal_id]['Thresholds'][index_spectrum] * self.crystals_maps[crystal_id]['Gain'][index_spectrum] + self.crystals_maps[crystal_id]['Offset'][index_spectrum], 0, 1200, color="green", label="Threshold", linestyle='--')
                    else:
                        axes[i_row+1][i_col+1].set_xlabel("ADU")
                        if self.if_show_threshold_in_spectrum and 'Threshold' in self.crystals_maps[crystal_id].keys():
                            axes[i_row+1][i_col+1].vlines(self.crystals_maps[crystal_id]['Thresholds'][index_spectrum], 0, 1200, linestyle='--', color="green", label="Threshold")
                    axes[i_row+1][i_col+1].set_ylabel("Counts")
                    axes[i_row+1][i_col+1].set_ylim(0, self.spectrum_ymax)
                    axes[i_row+1][i_col+1].legend()
        plt.show()

    def view_pixels_calibration(self, row, col, crystal_id):
        if "Gain" not in self.crystals_maps[crystal_id].keys() or "Offset" not in self.crystals_maps[crystal_id].keys():
            print(f"Calibration not found for crystal {crystal_id}")
            return
        fig, axes = plt.subplots(3, 3, figsize=(9, 9))
        plt.suptitle(f"Crystal {crystal_id} - Calibration")
        for i_row in [-1, 0, 1]:
            for i_col in [-1, 0, 1]:
                if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                    index_spectrum = (row + i_row) * self.n_col + (col + i_col)
                    energies = self.bins * self.crystals_maps[crystal_id]['Gain'][index_spectrum] + self.crystals_maps[crystal_id]['Offset'][index_spectrum]
                    axes[i_row+1][i_col+1].plot(self.bins, energies, '--')
                    for name, peaks in self.all_multiple_peaks[crystal_id].items():
                        if name != 'nosource':
                            axes[i_row+1][i_col+1].plot(peaks[index_spectrum], self.all_multiple_energies[crystal_id][name], color="red", linestyle='none', label=name, marker='+')
                    axes[i_row+1][i_col+1].set_title(f"Pixel ({row + i_row}, {col + i_col})")
                    axes[i_row+1][i_col+1].set_xlabel("ADU")
                    axes[i_row+1][i_col+1].set_ylabel("energy (keV)")
        plt.show()

    def view_pixels_scatters(self, row, col, crystal_id):
        fig, axes = plt.subplots(3, 3, figsize=(9, 9))
        plt.suptitle(f"Crystal {crystal_id} - Scatters")
        for i_row in [-1, 0, 1]:    
            for i_col in [-1, 0, 1]:

                if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                    
                    # load the scatter
                    pixel_id = (row + i_row) * self.n_col + (col + i_col)
                    print(f"Loading scatter of Pixel ID={pixel_id} ({row + i_row}, {col + i_col})")
                    scatters = self.lazy_load_scatters(row + i_row, col + i_col, crystal_id, axes[i_row+1][i_col+1])
                    if scatters is None:
                        continue

                    # draw the scatter
                    axes[i_row+1][i_col+1].plot(scatters, markersize=0.3, linestyle='none', marker='o')
                    axes[i_row+1][i_col+1].set_title(f"Pixel ({row + i_row}, {col + i_col})")
                    axes[i_row+1][i_col+1].set_xlabel("Frame")
                    axes[i_row+1][i_col+1].set_ylabel("ADU")
                    axes[i_row+1][i_col+1].legend()

                    # additional lines
                    # if self.if_show_peaks_in_scatters:
                    #     for name, peaks in self.all_multiple_peaks[crystal_id].items():
                    #         axes[i_row+1][i_col+1].hlines(peaks[row + i_row, col + i_col], 0, len(scatters), color="red", label=name)
                    if self.if_show_threshold_in_scatters:
                        axes[i_row+1][i_col+1].hlines(self.crystals_maps[crystal_id]['Threshold'][pixel_id], 0, len(scatters), color="green", label="Threshold")
        plt.show()

    def lazy_load_scatters(self, row, col, crystal_id, axes):
        if self.title == "CZT":
            # in matlab file convention, row and col start from 1, so we need to add 1
            scatter_file = f"{self.scatter_folder}\\scatter_{row + 1}_{col + 1}.mat"
            scatters = sio.loadmat(scatter_file)['scat'][0]
            return scatters

        if "alpha" in self.title.lower():
            pixel_id = row * self.n_col + col
            scatter_file = rf"{self.scatter_folder}\scatter_crystal_{crystal_id}_pixel_{pixel_id}.h5"

            if pathlib.Path(scatter_file).exists():
                with h5py.File(scatter_file, 'r') as f:
                    scatters = np.asarray(f['scatter'][:])
                    return scatters
            else:
                print(f"Scatter file not found: {scatter_file}")


