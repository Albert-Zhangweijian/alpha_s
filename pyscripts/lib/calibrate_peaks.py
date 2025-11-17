import h5py
import numpy as np
import matplotlib.pyplot as plt

def calibrate_peaks(peaks1, peaks2, energy1, energy2, threshold_energy=5, peaks1_start_adus=None, peaks2_start_adus=None, verbose=False):
    """
    Calibrate the peaks using linear fit (y = slope * x + intercept)
    peaks1, peaks2: arrays of peak positions in channel numbers
    energy1, energy2: values of energy corresponding to the peaks
    Returns: slope, intercept, R2, thresholds
    """

    slopes = (energy2 - energy1) / (peaks2 - peaks1)
    intercepts = energy1 - slopes * peaks1
    thresholds = np.floor((threshold_energy - intercepts) / slopes).astype(np.int32)
    R2 = np.ones_like(slopes)  # Placeholder for R2 values, can be computed if needed

    if (peaks1_start_adus is not None) and (peaks2_start_adus is not None):
        peaks1_start_erg = slopes * peaks1_start_adus + intercepts
        peaks2_start_erg = slopes * peaks2_start_adus + intercepts
        isvalids = (peaks1_start_erg < 10) & (peaks1_start_erg > -6) & (peaks2_start_erg < 10) & (peaks2_start_erg > -6)
    else:
        isvalids = np.ones_like(slopes, dtype=bool)

    return slopes, intercepts, isvalids, thresholds, R2

if __name__ == "__main__":

    for crystal_id in range(1, 4):
        peak1_filepath = rf"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_10hrs\panel_1\pixels_peaks_crystal_{crystal_id}.h5"
        peak2_filepath = rf"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_4hours\panel_1\pixels_peaks_crystal_{crystal_id}.h5"
        calibration_output_filepath = rf"H:\alpha_spect_mini\calibration_results\panel_1\calibration_crystal_{crystal_id}_new.h5"

        with h5py.File(peak1_filepath, 'r') as f:
            peaks1 = np.asarray(f['peak0'][:, 0])
            peaks1_start_adus = np.asarray(f['pixels_startadus'][:])
        with h5py.File(peak2_filepath, 'r') as f:
            peaks2 = np.asarray(f['peak0'][:, 0])
            peaks2_start_adus = np.asarray(f['pixels_startadus'][:])

        energy1 = 122.06  # keV for Co57
        energy2 = 59.54   # keV for Am241

        slope, intercept, isvalids, threshold, R2 = calibrate_peaks(peaks1, peaks2, energy1, energy2, \
                threshold_energy=5, peaks1_start_adus=peaks1_start_adus, peaks2_start_adus=peaks2_start_adus, verbose=True)
        
        plt.figure(figsize=(10, 10))
        plt.subplot(2,2,1)
        plt.imshow((peaks1_start_adus*slope+intercept).reshape(80,80), cmap='gray')
        plt.colorbar(label='Energy (keV)')
        plt.clim(-10, 10)
        plt.subplot(2,2,2)
        plt.imshow((peaks2_start_adus*slope+intercept).reshape(80,80), cmap='gray')
        plt.colorbar(label='Energy (keV)')
        plt.clim(-10, 10)
        plt.subplot(2,2,3)
        plt.imshow(isvalids.reshape(80,80), cmap='gray')
        plt.title(f"Panel 1, Crystal {crystal_id} - Valid Pixels Map")
        plt.colorbar(label='Is Valid Pixel')
        plt.show()
        
        with h5py.File(calibration_output_filepath, 'w') as f:
            f.create_dataset('pixels_calibrations', data=np.vstack((slope, intercept)).T)
            f.create_dataset('pixels_peaks', data=np.vstack((peaks1, peaks2)).T)
            f.create_dataset('pixels_isvalids', data=isvalids)
            f.create_dataset('pixels_thresholds', data=threshold)
