#include <array>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cstdint>
#include <stdexcept>
#include <cstdlib>
#include "highfive/HighFive.hpp"

#include "header.hpp"

int raw2frames(const std::vector<std::string> rawfiles, const std::string crystals_frames_folder) {

    // create output folder
    if (!std::filesystem::exists(crystals_frames_folder)) {
        std::filesystem::create_directories(crystals_frames_folder);
    }

    // define the data container for frames images
    uint16_t* crystals_frames_images = (uint16_t*) calloc(
        static_cast<size_t>(global_config["N_CRYSTALS"][0]) * global_config["MAX_FRAMES_SIZE"][0] * global_config["N_PIXELS"][0],
        sizeof(uint16_t)
    );

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
        std::cout << "Reading file: " << rawfile << ", total frames: " << n_frames << ", divided into " << n_flow_times << " flows." << std::endl;
        for (int flow_i = 0; flow_i < n_flow_times; flow_i++) {

            int start_frame_id = flow_i * global_config["MAX_FRAMES_SIZE"][0];
            int end_frame_id = start_frame_id + global_config["MAX_FRAMES_SIZE"][0] < n_frames ? start_frame_id + global_config["MAX_FRAMES_SIZE"][0] : n_frames;

            read_rawfile_to_crystals_frames_images(rawfile, start_frame_id, end_frame_id, crystals_frames_images);

            // save the frames using HighFive (HDF5 C++ wrapper)
            const size_t frames_in_chunk = static_cast<size_t>(end_frame_id - start_frame_id);
            const size_t n_pixels = static_cast<size_t>(global_config["N_PIXELS"][0]);
            for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; ++crystal_id) {
                std::string out_path = crystals_frames_folder + "\\" + format_string("frames_{}_{}_crystal_{}.h5", start_frame_id, end_frame_id, crystal_id);
                HighFive::File file(out_path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
                
                
                std::vector<size_t> dims{ static_cast<size_t>(global_config["MAX_FRAMES_SIZE"][0]), n_pixels };
                HighFive::DataSpace space(dims);
                auto dset = file.createDataSet<uint16_t>("frames", space);
                const uint16_t* data_ptr = crystals_frames_images + (static_cast<size_t>(crystal_id) * global_config["MAX_FRAMES_SIZE"][0]) * global_config["N_PIXELS"][0];
                dset.write_raw(data_ptr);
                // std::vector<size_t> offset = {0, 0};
                // std::vector<size_t> count  = {frames_in_chunk, n_pixels};
                // dset.select(offset, count).write(data_ptr);
            }

        }
    }



    free(crystals_frames_images);
    return 0;
}
