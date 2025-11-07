#include <vector>
#include <string>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#include "header.hpp"

int raw2scatters(const std::vector<std::string>& rawfiles, const std::string crystals_scatters_folder,
    const std::vector<std::string>& crystals_calibration_files, const std::string process_type) {


    // read calibration file, here we mainly use the pixels_isvalids to select pixels to check
    // besides, it also servers to mark the peak positions on the scatter plot.
    std::vector<std::vector<std::vector<float>>> crystals_pixels_peaks(global_config["N_CRYSTALS"][0], std::vector<std::vector<float>>(global_config["N_PIXELS"][0]));
    std::vector<std::vector<std::vector<float>>> crystals_pixels_calibration(global_config["N_CRYSTALS"][0], std::vector<std::vector<float>>(global_config["N_PIXELS"][0]));
    std::vector<std::vector<uint8_t>> crystals_pixels_isvalids(global_config["N_CRYSTALS"][0], std::vector<uint8_t>(global_config["N_PIXELS"][0], false));
    std::vector<std::vector<float>> crystals_pixels_thresholds(global_config["N_CRYSTALS"][0], std::vector<float>(global_config["N_PIXELS"][0]));
    for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
        if (crystals_calibration_files[crystal_id] == "skip") continue;
        read_pixels_calibrations(crystals_calibration_files[crystal_id], 
            crystals_pixels_calibration[crystal_id],
            crystals_pixels_peaks[crystal_id],
            crystals_pixels_isvalids[crystal_id],
            crystals_pixels_thresholds[crystal_id]);
    }

    // select pixels to check
    std::vector<std::vector<int>> crystals_target_ids(global_config["N_CRYSTALS"][0]);
    for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
        for (int pixel_id = 0; pixel_id < global_config["N_PIXELS"][0]; pixel_id++) {
            if (process_type == "all")
                crystals_target_ids[crystal_id].push_back(pixel_id);
            else if (process_type == "all-random" && (rand() % 100 < 5)) // if type is all-random, then randomly select 5% of pixels
                crystals_target_ids[crystal_id].push_back(pixel_id);
            else if (process_type == "valid-random" && crystals_pixels_isvalids[crystal_id][pixel_id] && (rand() % 100 < 5)) // randomly select 5% of the valid pixels
                crystals_target_ids[crystal_id].push_back(pixel_id);
            else if (process_type == "invalid" && !crystals_pixels_isvalids[crystal_id][pixel_id])
                crystals_target_ids[crystal_id].push_back(pixel_id);
        }
    }


    // create the folder and open each file
    if (!std::filesystem::exists(crystals_scatters_folder))
        std::filesystem::create_directory(crystals_scatters_folder);
    std::vector<std::vector<std::string>> crystals_pixels_scatter_files(global_config["N_CRYSTALS"][0]);
    for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
        for (int pixel_i = 0; pixel_i < crystals_target_ids[crystal_id].size(); pixel_i++) {
            int pixel_id = crystals_target_ids[crystal_id][pixel_i];
            std::string scatter_file = format_string("{}\\scatter_crystal_{}_pixel_{}.bin", crystals_scatters_folder, crystal_id, pixel_id);
            if (std::filesystem::exists(scatter_file)) {
                // std::cout << "Warning: scatter file already exists and will be overwritten: " << scatter_file << std::endl;
                std::filesystem::remove(scatter_file);
            }
            crystals_pixels_scatter_files[crystal_id].push_back(scatter_file);
        }
    }

    // define the data container for frames images (flat)
    const int n_crystals = global_config["N_CRYSTALS"][0];
    const int max_frames = global_config["MAX_FRAMES_SIZE"][0];
    const int n_pixels = global_config["N_PIXELS"][0];
    std::vector<uint16_t> crystals_frames_images(static_cast<size_t>(n_crystals) * max_frames * n_pixels, 0);

    // furthermore, since collapsing into scatters requires continuous loading of rawfiles into frames
    // the bottleneck is the reading speed of rawfiles which is already parallelized.
    // the process of converting frames into scatters is actually easy and not suitable for parallelization.
    // (expense for each pixel is small, but we can distribute multiple pixels into each thread) -- to be implemented in the future.

    const size_t n_frame_pixels = global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];
    for (auto& rawfile : rawfiles) {

        if (!std::filesystem::exists(rawfile))
            std::cout << "Warning: rawfile does not exist and will be skipped: " << rawfile << std::endl;
        if (!std::filesystem::is_regular_file(rawfile))
            throw std::runtime_error("ERROR: Not a regular file: " + rawfile);
        std::ifstream file(rawfile, std::ios::binary | std::ios::ate);
        if (!file.is_open()) throw std::runtime_error("Cannot open: " + rawfile);
        const size_t file_size = file.tellg();
        const size_t n_frames = file_size / n_frame_pixels / sizeof(uint16_t);
        file.close();

        int n_flow_times = (n_frames + global_config["MAX_FRAMES_SIZE"][0] - 1) / global_config["MAX_FRAMES_SIZE"][0];
        std::cout << "Reading file: " << rawfile << ", total frames: " << n_frames << ", divided into " << n_flow_times << " flows." << std::endl;
        for (int flow_i = 0; flow_i < n_flow_times; flow_i++) {
            int start_frame_id = flow_i * global_config["MAX_FRAMES_SIZE"][0];
            int end_frame_id = start_frame_id + global_config["MAX_FRAMES_SIZE"][0] < n_frames ? start_frame_id + global_config["MAX_FRAMES_SIZE"][0] : n_frames;
            read_rawfile_to_crystals_frames_images(rawfile, start_frame_id, end_frame_id, crystals_frames_images.data());

            const int frames_in_chunk = end_frame_id - start_frame_id;
            for (int crystal_id = 0; crystal_id < n_crystals; crystal_id++) {
                for (int pixel_i = 0; pixel_i < crystals_target_ids[crystal_id].size(); pixel_i++) {
                    int pixel_id = crystals_target_ids[crystal_id][pixel_i];

                    std::string scatter_file = crystals_pixels_scatter_files[crystal_id][pixel_i];
                    std::ofstream scatter_ostream(scatter_file, std::ios::binary | std::ios::app);
                    const uint16_t* frames_ptr = crystals_frames_images.data() + (static_cast<size_t>(crystal_id) * max_frames) * n_pixels;
                    for (int frame_i = 0; frame_i < frames_in_chunk; frame_i++) {
                        const uint16_t* frame_ptr = frames_ptr + static_cast<size_t>(frame_i) * n_pixels;
                        uint16_t value = frame_ptr[pixel_id];
                        if (value >= global_config["ADU_MAX"][0] || value < global_config["ADU_MIN"][0]) continue;
                        scatter_ostream.write(reinterpret_cast<const char*>(&value), sizeof(uint16_t));
                    }
                }
            }
        }
    }


    return 0;
}
