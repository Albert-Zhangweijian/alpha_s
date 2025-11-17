import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from scipy.special import erf
from scipy.optimize import approx_fprime

import warnings
warnings.filterwarnings("ignore", message="delta_grad == 0.0.*")

def gaussian(x, pos, height, sigma):
    return height * np.exp(-0.5 * ((x - pos) / sigma) ** 2)

def asymmetric_gaussian(x, pos, height, left_sigma, right_sigma, left_tail_sigma, dividing_pos):
    y = np.zeros_like(x)
    left_mask = x < pos
    right_mask = x >= pos
    tail_mask = (x < pos - dividing_pos * left_sigma) & left_mask
    tail_height = height * np.exp(-dividing_pos * dividing_pos / 2) * np.exp(dividing_pos * left_sigma * left_sigma * dividing_pos / 2 / left_tail_sigma / left_tail_sigma)
    y[left_mask] = height * np.exp(-0.5 * ((x[left_mask] - pos) / left_sigma) ** 2)
    y[tail_mask] = tail_height * np.exp(-0.5 * ((x[tail_mask] - pos) / left_tail_sigma) ** 2)
    y[right_mask] = height * np.exp(-0.5 * ((x[right_mask] - pos) / right_sigma) ** 2)
    return y

def convolved_erf_gaussian(x, pos, height, left_sigma, B, Gamma):
    y = height * np.exp(-0.5 * ((x - pos) / left_sigma) ** 2) + (1 - erf( (x - pos) / (np.sqrt(2) * left_sigma) + left_sigma / Gamma / np.sqrt(2) )) * B * np.exp((x - pos) / Gamma)
    return y

def double_convolved_erf_gaussian(x, pos, height, left_sigma, B, Gamma, B2, Gamma2):
    y = height * np.exp(-0.5 * ((x - pos) / left_sigma) ** 2) + \
        (1 - erf( (x - pos) / (np.sqrt(2) * left_sigma) + left_sigma / Gamma / np.sqrt(2) )) * B * np.exp((x - pos) / Gamma) + \
            (1 - erf( (x - pos) / (np.sqrt(2) * left_sigma) + left_sigma / Gamma2 / np.sqrt(2) )) * B2 * np.exp((x - pos) / Gamma2)
    return y

def generate_initial_guess_and_bounds(bins, spectrum, model):
    gaussian_params = [bins[np.argmax(spectrum)]+4, np.max(spectrum), (bins[-1] - bins[0]) / 8]  # pos, height, sigma
    if model == 'asymmetric_gaussian':
        initial_guess = [gaussian_params[0], gaussian_params[1], gaussian_params[2]*1.5, gaussian_params[2], gaussian_params[2]*4, 1.25]
        upper_bounds = [gaussian_params[0] + 16, gaussian_params[1] * 2, 100, 100, 100, 2.0]
        lower_bounds = [gaussian_params[0] - 16, gaussian_params[1] / 2, 7, 10, 10, 1.0]
        constraints = [
            {'type': 'ineq', 'fun': lambda x:  x[2] - 1.05*x[3]},
            {'type': 'ineq', 'fun': lambda x:  x[4] - 1.25*x[2]}
        ]
        penalty = lambda params: 0.01 * (2*params[2] - params[4])**2
    elif model == 'convolved_erf_gaussian':
        initial_guess = [gaussian_params[0], gaussian_params[1]*3/4, gaussian_params[2] / 2, gaussian_params[1]/4, 80]
        upper_bounds = [gaussian_params[0] + 36, gaussian_params[1], 20, gaussian_params[1], 100]
        lower_bounds = [gaussian_params[0] - 36, gaussian_params[1] / 5, 1, 1, 10]
        constraints = []
        penalty = lambda params: 0
    elif model == 'double_convolved_erf_gaussian':
        initial_guess = [gaussian_params[0], gaussian_params[1]*3/4, gaussian_params[2] / 2, gaussian_params[1]/4, 40, gaussian_params[1]/8, 150]
        upper_bounds = [gaussian_params[0] + 36, gaussian_params[1], 30, gaussian_params[1]/2, 100, gaussian_params[1]/2, 300]
        lower_bounds = [gaussian_params[0] - 36, gaussian_params[1]/5, 1, 1, 10, 1, 40]
        constraints = []
        penalty = lambda params: 0
    else:
        raise ValueError("Unsupported model type: {}".format(model))
    initial_guess = np.array(initial_guess)
    lower_bounds = np.array(lower_bounds)
    upper_bounds = np.array(upper_bounds)
    return initial_guess, lower_bounds, upper_bounds, constraints, penalty


def fit_peaks(bins, spectrum, peak_index, fit_range_left, fit_range_right, smooth=True, verbose=False, model='asymmetric_gaussian'):
    
    if model == 'asymmetric_gaussian':
        model_function = asymmetric_gaussian
    elif model == 'convolved_erf_gaussian':
        model_function = convolved_erf_gaussian
    elif model == 'double_convolved_erf_gaussian':
        model_function = double_convolved_erf_gaussian
    else:
        raise ValueError("Unsupported model type: {}".format(model))

    left_boundary = max(0, peak_index - fit_range_left)
    right_boundary = min(len(bins) - 1, peak_index + fit_range_right)
    spectrum = np.convolve(spectrum, [1/4, 1/2, 1/4], mode='same') if smooth else spectrum
    partial_bins = bins[left_boundary:right_boundary+1]
    partial_spectrum = spectrum[left_boundary:right_boundary+1]

    initial_guess, lower_bounds, upper_bounds, constraints, penalty = generate_initial_guess_and_bounds(partial_bins, partial_spectrum, model=model)
    objective = lambda params: np.sum((partial_spectrum - model_function(partial_bins, *(params)))**2) + penalty(params)
    bounds = opt.Bounds(lower_bounds, upper_bounds)

    
    normalized_initial_guess = (initial_guess - lower_bounds) / (upper_bounds - lower_bounds) * 1000
    objective = lambda params: np.sum((partial_spectrum - model_function(partial_bins, *(params*(upper_bounds - lower_bounds)  / 1000 + lower_bounds)))**2) + penalty(params*(upper_bounds - lower_bounds) / 1000 + lower_bounds)
    bounds = opt.Bounds([1e-5]*len(lower_bounds), [1000-1e-5]*len(upper_bounds))

    result = opt.minimize(objective, normalized_initial_guess, bounds=bounds, method='COBYLA', options={'maxiter': 5000, 'disp': False, 'xtol': 1e-12, 'gtol': 1e-8}, constraints=constraints)
    # result = opt.minimize(objective, normalized_initial_guess, bounds=bounds, method='COBYLA', options={'maxiter': 5000, 'disp': False, 'xtol': 1e-22, 'gtol': 1e-12}, constraints=constraints)
    fitted_params = result.x / 1000 * (upper_bounds - lower_bounds) + lower_bounds

    # In principle, there should always be a result
    # But the correctness of the result needs to be verified after the fitting
    if not result.success:
        print(f"Warning: Fitting did not converge for peak at index {peak_index}. Message: {result.message}. But proceeding with the best available fit.")

    if verbose:

        # Print fitting results
        print(f"Fitting result with success={result.success}, loss={objective(result.x):.2f}, params={fitted_params}, message: {result.message}")
        print("Gradient at minimum:", approx_fprime(result.x, objective, 1e-8))
        if model == 'asymmetric_gaussian':
            print(f"Final/Intial: Left_sigma={fitted_params[2]/initial_guess[2]:.2f} Right_sigma={fitted_params[3]/initial_guess[3]:.2f} Left_tail_sigma={fitted_params[4]/initial_guess[4]:.2f}")
            print(f"Fitted Sigma Ratios:  Left_sigma/right_sigma={fitted_params[2]/fitted_params[3]:.2f} Left_tail_sigma/right_sigma={fitted_params[4]/fitted_params[3]:.2f}")
        
        partial_fine_bins = np.linspace(partial_bins[0], partial_bins[-1], 200)
        plt.figure(figsize=(8, 5)) 
        
        plt.plot(bins, spectrum, label='Original Data', marker='o', markersize=4, linestyle='None', alpha=0.5)
        plt.plot(bins[peak_index], spectrum[peak_index], marker='o', color='red', linestyle='None', label='Detected Peak')
        plt.axvline(partial_bins[0], color='green', linestyle='--', label='Fitting Range')
        plt.axvline(partial_bins[-1], color='green', linestyle='--')

        initial_curve = model_function(partial_fine_bins, *initial_guess)
        fitted_curve = model_function(partial_fine_bins, *fitted_params)
        plt.plot(partial_fine_bins, initial_curve, label='Initial Guess', linestyle='--', color='green', alpha=0.7)
        plt.plot(partial_fine_bins, fitted_curve, label='Fitted Curve', linestyle='-', color='green')
        if model == 'asymmetric_gaussian':
            plt.plot(fitted_params[0] - fitted_params[2] * fitted_params[5], fitted_params[1] * np.exp(-0.5 * (fitted_params[5] ** 2)), marker='x', linestyle='None', color='purple', label='Tail Dividing Point')
        if model in ['convolved_erf_gaussian', 'double_convolved_erf_gaussian']:
            plt.plot(partial_fine_bins, gaussian(partial_fine_bins, *fitted_params[:3]), label='Gaussian Component', linestyle='--', color='gray')
            plt.plot(partial_fine_bins, model_function(partial_fine_bins, *fitted_params) - gaussian(partial_fine_bins, *fitted_params[:3]), label='Tail Component', linestyle=':', color='gray')

        plt.title('Peak Fitting')
        plt.xlabel('Bins')
        plt.ylabel('Counts')
        plt.ylim(0, np.max(partial_spectrum)*1.2)
        plt.xlim(partial_bins[0]-100, partial_bins[-1]+100)
        plt.legend()
        plt.grid()
        # plt.show()

    return fitted_params, objective(result.x)

if __name__ == "__main__":

    import h5py
    
    import os, sys
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from lib.find_peaks import find_peaks

    n_peaks = 1
    spectra_file = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\crystals_spectra\pixels_spectra_crystal_0.h5"
    peaks_file = r"H:\alpha_spect_mini\20241120_EneCalib_Co57_150fps_-500V_7hrs\panel_1\pixels_peaks_crystal_0.h5"
    # spectra_file = r"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_4hours\panel_1\pixels_spectra_crystal_0.h5"
    # peaks_file = r"H:\alpha_spect_mini\20241220_Am241_panel12356_800fps_4hours\panel_1\pixels_peaks_crystal_0.h5"
    with h5py.File(spectra_file, 'r') as f:
        bins = f['bins'][:]
        pixels_spectra = f['pixels_spectra'][:]
    with h5py.File(peaks_file, 'r') as f:
        pixels_peaks = f['peak0'][:]


    for pixel_id in range(730, pixels_spectra.shape[0]):
        print(f"Processing pixel {pixel_id}/{pixels_spectra.shape[0]}")
        peaks_indices = find_peaks(bins, pixels_spectra[pixel_id:pixel_id+1, :], n_peaks=n_peaks, kernel_half_length=3, neglected_bins=200, peak_distance=5, verbose=False)

        fit_range_left = 6
        fit_range_right = 6
        for peak_i in range(1):
            peak_index = peaks_indices[peak_i]
            fitted_params = fit_peaks(bins, pixels_spectra[pixel_id], peak_index, fit_range_left, fit_range_right, smooth=False, verbose=True, model='asymmetric_gaussian')
            partial_fine_bins = np.linspace(bins[peak_index - fit_range_left], bins[peak_index + fit_range_right], 200)
            
            plt.gcf()

            # # smoothed fitted curve
            # smooth_fitted_params = fit_peaks(bins, pixels_spectra[pixel_id], peak_index, fit_range_left, fit_range_right, smooth=True, verbose=False)
            # plt.plot(partial_fine_bins, asymmetric_gaussian(partial_fine_bins, *smooth_fitted_params), label='Fitted Curve(Smoothed)', linestyle='-')

            # old fitted curve
            # fitted_curve = asymmetric_gaussian(partial_fine_bins, *pixels_peaks[pixel_id, :])
            # plt.plot(partial_fine_bins, fitted_curve, label='Fitted Curve(Old)', linestyle='--', color='orange')

            plt.title(f'Pixel {pixel_id}, Peak {peak_i}')
            plt.legend()
            plt.show()