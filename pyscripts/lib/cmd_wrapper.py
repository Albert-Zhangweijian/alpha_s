import subprocess

def run_cmd(args, kwargs):
    args_config = [] if kwargs.get('config') is None else ['--config', kwargs.get('config')]
    full_args = ["alpha.exe"] + args_config + args
    cwd = '.' if kwargs.get('cwd') is None else str(kwargs.get('cwd'))
    print("===== Running command: =====")
    print("Config:", args_config)
    print("Working directory:", cwd)
    print("Full command args:", full_args)
    print("============================")

    subprocess.run(full_args, cwd=cwd)
    # result = subprocess.run(args, cwd=cwd, stderr=subprocess.PIPE, text=True)
    # if result.returncode != 0:
    #     print(f"Error occurred: {result.stderr}")
    # return result

def raw2frames(rawfiles, frames_folder, **kwargs):
    args = ["raw2frames"] \
    + ["-i"] + rawfiles \
    + ["-o", frames_folder]
    run_cmd(args, kwargs)

def raw2spectra(rawfiles, crystals_spectra_folder, **kwargs):

    args = ["raw2spectra"] \
    + ["-i"] + rawfiles \
    + ["-o", crystals_spectra_folder]
    run_cmd(args, kwargs)

def raw2scatters(rawfiles, crystals_scatters_folder, crystals_calibrations, mode, **kwargs):
    args = ["raw2scatters"] \
        + ["-i"] + rawfiles \
        + ["-o", crystals_scatters_folder] \
        + [elem for calib in crystals_calibrations for elem in ("-c", calib)] \
        + ["-m", mode]
    run_cmd(args, kwargs)

def raw2clusters(rawfiles, crystals_clusters_folder, crystals_calibrations, crystals_thresholds, extended_mode, **kwargs):
    args = ["raw2clusters"] \
    + ["-i"] + rawfiles \
    + ["-o", crystals_clusters_folder] \
    + [elem for crystal_calibration in crystals_calibrations for elem in ("-c", crystal_calibration)] \
    + [elem for crystal_threshold in crystals_thresholds for elem in ("-t", crystal_threshold)] \
    + (["--extended_mode"] if extended_mode else [])
    run_cmd(args, kwargs)