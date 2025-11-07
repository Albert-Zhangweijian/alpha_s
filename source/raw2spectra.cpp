#include <array>
#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include "highfive/HighFive.hpp"

#include "header.hpp"

std::vector<size_t> compute_adu_lookup_table(const std::vector<float>& bins) {
    std::vector<size_t> adu_lookup_table(global_config["ADU_MAX"][0] - global_config["ADU_MIN"][0], 0);
    for (int adu = global_config["ADU_MIN"][0]; adu < global_config["ADU_MAX"][0]; ++adu) {
        auto it = std::upper_bound(bins.begin(), bins.end(), adu);
        adu_lookup_table[adu-global_config["ADU_MIN"][0]] = (it != bins.begin()) ? std::distance(bins.begin(), it) - 1 : 0;
    }
    return adu_lookup_table;
}


int read_frames_images_to_pixels_spectra(
    const uint16_t* frames_images,
    const std::vector<float>& /*bins*/,
    const std::vector<size_t>& lookup_table,
    int* pixels_spectra,
    int n_frames) {

    const int n_pixels = global_config["N_PIXELS"][0];
    const int adu_min = global_config["ADU_MIN"][0];
    const int adu_max = global_config["ADU_MAX"][0];
    const int n_bins = global_config["N_BINS"][0];

    #pragma omp parallel for
    for (int frame_i = 0; frame_i < n_frames; ++frame_i) {
        const uint16_t* frame_image = frames_images + static_cast<size_t>(frame_i) * n_pixels;
        for (int pixel_id = 0; pixel_id < n_pixels; ++pixel_id) {
            uint16_t value = frame_image[pixel_id];
            if (value >= adu_max || value < adu_min) continue;
            size_t bin_idx = lookup_table[static_cast<size_t>(value - adu_min)];
            if (bin_idx >= static_cast<size_t>(n_bins)) continue;
            pixels_spectra[static_cast<size_t>(pixel_id) * n_bins + bin_idx]++;
        }
    }
    return 0;
}

int raw2spectra(const std::vector<std::string> rawfiles, const std::string crystals_spectra_folder) {

    // create output folder
    if (!std::filesystem::exists(crystals_spectra_folder)) {
        std::filesystem::create_directories(crystals_spectra_folder);
    }

    // create data containers
    const int n_crystals = global_config["N_CRYSTALS"][0];
    const int n_pixels = global_config["N_PIXELS"][0];
    const int n_bins = global_config["N_BINS"][0];
    const int max_frames = global_config["MAX_FRAMES_SIZE"][0];
    std::vector<int> crystals_pixels_spectra(static_cast<size_t>(n_crystals) * n_pixels * n_bins, 0);

    // precompute the lookup table for adu to bin index
    std::vector<float> bins = create_bins(global_config["N_BINS"][0], global_config["ADU_MIN"][0], global_config["ADU_MAX"][0]);
    static const std::vector<size_t> lookup_table = compute_adu_lookup_table(bins);

    // define the data container for frames images (flat)
    std::vector<uint16_t> crystals_frames_images(static_cast<size_t>(n_crystals) * max_frames * n_pixels, 0);

    // collapse each raw data file into a spectrum
    // In order to limit the size of these frames, they should be updated for every file;
    // otherwise, they will be very large after reading all files.
    // Everytime a file is read, we tally it immediately.
    const size_t n_frame_pixels = global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];
    for (auto& rawfile : rawfiles) {

        // since sometimes we deal with ultra-large rawfiles, we must read and analyze them in chunks.
        std::ifstream file(rawfile, std::ios::binary | std::ios::ate);
        if (!file.is_open()) throw std::runtime_error("Cannot open: " + rawfile);
        const size_t file_size = file.tellg();
        const size_t n_frames = file_size / n_frame_pixels / sizeof(uint16_t);
        file.close();

        // divide 1 file into multiple chunks and then read in parallel
        int n_flow_times = (n_frames + global_config["MAX_FRAMES_SIZE"][0] - 1) / global_config["MAX_FRAMES_SIZE"][0];
        for (int flow_i = 0; flow_i < n_flow_times; flow_i++) {

            int start_frame_id = flow_i * global_config["MAX_FRAMES_SIZE"][0];
            int end_frame_id = start_frame_id + global_config["MAX_FRAMES_SIZE"][0] < n_frames ? start_frame_id + global_config["MAX_FRAMES_SIZE"][0] : n_frames;

            read_rawfile_to_crystals_frames_images(rawfile, start_frame_id, end_frame_id, crystals_frames_images.data());

            const int frames_in_chunk = end_frame_id - start_frame_id;
            for (int crystal_id = 0; crystal_id < n_crystals; ++crystal_id) {
                const uint16_t* frames_ptr = crystals_frames_images.data() + (static_cast<size_t>(crystal_id) * max_frames) * n_pixels;
                int* spectra_ptr = crystals_pixels_spectra.data() + (static_cast<size_t>(crystal_id) * n_pixels) * n_bins;
                read_frames_images_to_pixels_spectra(frames_ptr, bins, lookup_table, spectra_ptr, frames_in_chunk);
            }

        }
    }

    //save the spectrum
    for (int crystal_id = 0; crystal_id < n_crystals; ++crystal_id) {
        std::string out_path = crystals_spectra_folder + "\\" + format_string("pixels_spectra_crystal_{}.h5", crystal_id);
        HighFive::File file(out_path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        // bins 1D (simple API)
        file.createDataSet("bins", bins);

        // pixels_spectra 2D [N_PIXELS x N_BINS] from flat buffer
        std::vector<size_t> dims{ static_cast<size_t>(n_pixels), static_cast<size_t>(n_bins) };
        HighFive::DataSpace space(dims);
        auto dset = file.createDataSet<int>("pixels_spectra", space);
        const int* data_ptr = crystals_pixels_spectra.data() + (static_cast<size_t>(crystal_id) * n_pixels) * n_bins;
        dset.write_raw(data_ptr);
    }

    return 0;
}
