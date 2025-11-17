import h5py
import numpy as np
import matplotlib.pyplot as plt


def show_fitting_plot(event, bins, pixels_spectra, model, pixels_peak0_params, pixels_peak1_params):
    toolbar = plt.get_current_fig_manager().toolbar
    if toolbar.mode != '':  # when using zoom, pan, etc. ignore clicks
        print(f"Click ignored: {toolbar.mode} mode is active.")
        return

    if event.xdata is None or event.ydata is None:
        return

    row, col = round(event.ydata), round(event.xdata)
    print(f"Clicked on ({event.ydata:.2f}, {event.xdata:.2f}) => Row={row} Col={col}")
    pixel_peak0_params = pixels_peak0_params[row*40 + col]
    pixel_peak1_params = pixels_peak1_params[row*40 + col]
    plt.figure()
    plt.subplot(1, 3, 1)
    plt.plot(bins, pixels_spectra[row*40 + col], label='Spectrum', linestyle='none', marker='.', color='blue', markersize=4)
    plt.axvline(x=pixel_peak0_params[0], color='orange', linestyle='--', label='Peak 0 Mean')
    plt.axvline(x=pixel_peak1_params[0], color='green', linestyle='--', label='Peak 1 Mean')
    plt.xlim([0, 6000])
    plt.ylim([0, 500])
    plt.subplot(1, 3, 2)
    plt.plot(bins, pixels_spectra[row*40 + col], label='Spectrum', linestyle='none', marker='.', color='blue', markersize=4)
    plt.plot(bins, model(bins, *pixel_peak0_params[:6]), label='Fitted Peak 0', linestyle='-', color='orange')
    plt.xlim(0, pixel_peak0_params[0] + 200)
    plt.ylim(0, pixel_peak0_params[1]*1.5)
    plt.title(f'Loss = {pixel_peak0_params[6]:.2f}')
    plt.subplot(1, 3, 3)
    plt.plot(bins, pixels_spectra[row*40 + col], label='Spectrum', linestyle='none', marker='.', color='blue', markersize=4)
    plt.plot(bins, model(bins, *pixel_peak1_params[:6]), label='Fitted Peak 1', linestyle='-', color='green')
    plt.xlim(0, pixel_peak1_params[0] + 200)
    plt.ylim(0, pixel_peak1_params[1]*1.5)
    plt.title(f'Loss = {pixel_peak1_params[6]:.2f}')
    plt.show()

def interactive_calibration_map(bins, pixels_spectra, model, pixels_peak0_params, pixels_peak1_params):

    fig = plt.figure(figsize=(8, 6))
    ax = plt.imshow(pixels_peak0_params[:,1].reshape(40,40), cmap='gray')
    plt.colorbar(label='Peak 0 Mean (ADU)')
    plt.title('Peak 0 Mean Map')
    fig.canvas.mpl_connect('button_press_event', lambda event: show_fitting_plot(event, bins, pixels_spectra, model, pixels_peak0_params, pixels_peak1_params))
    plt.show()
