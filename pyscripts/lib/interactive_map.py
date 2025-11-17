import h5py
import numpy as np
import scipy.io as sio
import pathlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection



class InteractiveMap:

    def __init__(self, map, pixels_spectra, maprange=(0, 1), shape=(80, 80), adu_bins=None, calibration_params=None):
        self.map = map
        self.pixels_spectra = pixels_spectra
        self.shape = shape
        self.n_row = shape[0]
        self.n_col = shape[1]
        self.range = maprange
        self.spectrum_range = None

        if adu_bins is not None:
            self.adu_bins = adu_bins
        else:
            self.adu_bins = np.linspace(6000, 14000, 1001)[:-1]
        if calibration_params is not None:
            self.calibration_params = calibration_params
        else:
            self.calibration_params = np.zeros((self.n_row*self.n_col, 2))
            self.calibration_params[:, 0] = 1.0

        # draw map
        self.fig = plt.figure()
        self.ax = plt.gca()
        im = self.ax.imshow(self.map.reshape(self.shape), cmap="hot", 
                            origin="upper", interpolation='none',
                            vmin=self.range[0], vmax=self.range[1])
        plt.colorbar(im, ax=self.ax)
        self.ax.set_aspect('equal')
        self.ax.set_xlabel("Pixel X")
        self.ax.set_ylabel("Pixel Y")

        # add gridlines
        grid_lines = []
        for i in range(self.n_row):
            grid_lines.append([(-0.5, i+0.5), (self.n_col-0.5, i+0.5)])
        for j in range(self.n_col):
            grid_lines.append([(j+0.5, -0.5), (j+0.5, self.n_row-0.5)])
        grid = LineCollection(grid_lines, color="black", linewidths=0.5)
        self.ax.add_collection(grid)

        # add click event
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
        self.ax.plot(col, row, 'cs', fillstyle='none', markersize=4)[0]
        print(f"Clicked on ({event.ydata:.2f}, {event.xdata:.2f}) => Row={row} Col={col}")
        self.view_pixels_spectrum(row, col)

    def view_pixels_spectrum(self, row, col):
        fig, axes = plt.subplots(3, 3, figsize=(9, 9))
        
        for i_row in [-1, 0, 1]:
            for i_col in [-1, 0, 1]:
                if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                    pixel_id = (row + i_row) * self.n_col + (col + i_col)
                    xlabel = "Energy (keV)" if self.calibration_params[pixel_id, 1] != 0 else "ADU"
                    bins = self.adu_bins * self.calibration_params[pixel_id, 0] + self.calibration_params[pixel_id, 1]
                    axes[i_row+1][i_col+1].plot(bins, self.pixels_spectra[pixel_id, :], label='Spectrum')
                    axes[i_row+1][i_col+1].set_title(f"Pixel ({row + i_row}, {col + i_col}) ")
                    axes[i_row+1][i_col+1].set_ylabel("Counts")  
                    axes[i_row+1][i_col+1].set_xlabel(xlabel)
                    if self.spectrum_range is not None:
                        axes[i_row+1][i_col+1].set_xlim(self.spectrum_range)
                    axes[i_row+1][i_col+1].set_ylim(0, 10)
                    axes[i_row+1][i_col+1].legend()
        self.view_pixels_spectrum_callback(row, col, fig, axes)
        plt.show()

    def view_pixels_spectrum_callback(self, row, col, fig, axes):
        pass


class ThresholdMap(InteractiveMap):
    
    def __init__(self, *args, thresholds, **kargs):
        super().__init__(*args, **kargs)
        self.thresholds = thresholds

    def view_pixels_spectrum_callback(self, row, col, fig, axes):
        for i_row in [-1, 0, 1]:
            for i_col in [-1, 0, 1]:
                if (row + i_row) >= 0 and (row + i_row) < self.n_row and (col + i_col) >= 0 and (col + i_col) < self.n_col:
                    pixel_id = (row + i_row) * self.n_col + (col + i_col)
                    threshold = self.thresholds[pixel_id] * self.calibration_params[pixel_id, 0] + self.calibration_params[pixel_id, 1]
                    hrange = 100 * self.calibration_params[pixel_id, 0]
                    axes[i_row+1][i_col+1].vlines(threshold, 0, 30000, color='orange', linestyle='--')
                    # axes[i_row+1][i_col+1].set_xlim(threshold-hrange, threshold+hrange)
                    axes[i_row+1][i_col+1].set_ylim(0, 30000)
                    # axes[i_row+1][i_col+1].set_yscale('log')