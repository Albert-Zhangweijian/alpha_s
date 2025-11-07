#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include "header.hpp"

struct Cluster {
    int frame_id;
    std::vector<int> pixel_ids;
    std::vector<uint16_t> adus;
    std::vector<float> energies;
    // bool selected;
};


// For one pixel, add its adjacent pixels to the list of adjacent pixels if they are not already checked
int add_adjacent_pixels(int pixel_id, std::vector<int> &adjacent_pixel_ids, const std::vector<bool> &pixels_ischeckeds, bool include_diagonal) {

    std::vector<int> candidate_adjacent_pixel_ids;
    if (include_diagonal) {
        if (pixel_id % global_config["N_COLS"][0] == 0) {
            candidate_adjacent_pixel_ids = { pixel_id + 1, 
                pixel_id + global_config["N_COLS"][0], pixel_id + global_config["N_COLS"][0] + 1, 
                pixel_id - global_config["N_COLS"][0], pixel_id - global_config["N_COLS"][0] + 1 };
        } else if (pixel_id % global_config["N_COLS"][0] == global_config["N_COLS"][0] - 1) {
            candidate_adjacent_pixel_ids = { pixel_id - 1,
                pixel_id + global_config["N_COLS"][0], pixel_id + global_config["N_COLS"][0] - 1, 
                pixel_id - global_config["N_COLS"][0], pixel_id - global_config["N_COLS"][0] - 1 };
        } else 
            candidate_adjacent_pixel_ids = { pixel_id - 1, pixel_id + 1, 
                pixel_id + global_config["N_COLS"][0], pixel_id + global_config["N_COLS"][0] - 1, pixel_id + global_config["N_COLS"][0] + 1, 
                pixel_id - global_config["N_COLS"][0], pixel_id - global_config["N_COLS"][0] - 1, pixel_id - global_config["N_COLS"][0] + 1 };
    } else {
        if (pixel_id % global_config["N_COLS"][0] == 0) {  // pixel at the left edge, next pixel is the right edge of previous row
            candidate_adjacent_pixel_ids = { pixel_id + 1, 
                pixel_id + global_config["N_COLS"][0],
                pixel_id - global_config["N_COLS"][0] };
        } else if (pixel_id % global_config["N_COLS"][0] == global_config["N_COLS"][0] - 1) {  // pixel at the right edge, next pixel is the left edge of next row
            candidate_adjacent_pixel_ids = { pixel_id - 1,
                 pixel_id - global_config["N_COLS"][0], 
                 pixel_id + global_config["N_COLS"][0] };
        } else
            candidate_adjacent_pixel_ids = { pixel_id - 1, pixel_id + 1, 
                pixel_id + global_config["N_COLS"][0],
                pixel_id - global_config["N_COLS"][0] };
    }

    for (auto new_adjacent_pixel_id : candidate_adjacent_pixel_ids) {
        if (new_adjacent_pixel_id >= 0 && new_adjacent_pixel_id < global_config["N_PIXELS"][0] && !pixels_ischeckeds[new_adjacent_pixel_id])
            adjacent_pixel_ids.push_back(new_adjacent_pixel_id);
    }

    return 0;
}


std::vector<Cluster> read_frames_images_to_clusters(const uint16_t* frames_images,
    int n_frames,
    const float* calib0,
    const float* calib1,
    const uint8_t* pixels_isvalids,
    const float* pixels_thresholds,
    const float* pixels_secondary_thresholds,
    const int& global_start_frame_id,
    const bool extended_mode) {

    // distribute the frames to different threads
    int n_frames_per_thread = n_frames / global_config["N_THREADS"][0];
    std::vector<int> start_frame_ids(global_config["N_THREADS"][0]);
    for (int thread_i = 0; thread_i < global_config["N_THREADS"][0]; thread_i++)
        start_frame_ids[thread_i] = thread_i * n_frames_per_thread;
    start_frame_ids.push_back(n_frames);

    // find clusters in each frame in parallel
    std::vector<std::vector<Cluster>> threads_clusters(global_config["N_THREADS"][0]);
    #pragma omp parallel num_threads (global_config["N_THREADS"][0]) shared(start_frame_ids, calib0, calib1, pixels_isvalids, pixels_thresholds)
    {
        auto start_time = std::chrono::high_resolution_clock::now();

        // get the start and end frame ids for this thread
        int thread_id = omp_get_thread_num();
        int start_frame_id = start_frame_ids[thread_id];
        int end_frame_id = start_frame_ids[thread_id + 1];

        std::vector<bool> pixels_ischeckeds(global_config["N_PIXELS"][0], false);
        for (int frame_id = start_frame_id; frame_id < end_frame_id; frame_id++) {

            // get the frame image and set the bad pixels to be already checked
            const uint16_t* frame_image = frames_images + static_cast<size_t>(frame_id) * global_config["N_PIXELS"][0];
            for (int pixel_id = 0; pixel_id < global_config["N_PIXELS"][0]; pixel_id++)
                pixels_ischeckeds[pixel_id] = !pixels_isvalids[pixel_id] || (frame_image[pixel_id] < pixels_thresholds[pixel_id]);

            // check each pixel in order
            for (int pixel_id = 0; pixel_id < global_config["N_PIXELS"][0]; pixel_id++) {

                // though this pixel may be lower than the primary threshold
                // it may be included in the cluster because it can be larger than the secondary threshold
                // so we do not set it to be checked
                if (pixels_ischeckeds[pixel_id])
                    continue;

                // std::cout << pixel_id << " " << frame_image[pixel_id] << " " << pixels_thresholds[pixel_id] << " " << pixels_secondary_thresholds[pixel_id] << std::endl;

                // we first locate the center pixel, then put it into adjacent_pixel_ids
                // when the vector is not empty, we will pop the first element and check it
                // during checking, all other possible pixels will be added to the end of the vector
                std::vector<int> active_pixel_ids = {};
                std::vector<int> adjacent_pixel_ids = { pixel_id };
                while (!adjacent_pixel_ids.empty()) {
                    int adjacent_pixel_id = adjacent_pixel_ids.front();
                    adjacent_pixel_ids.erase(adjacent_pixel_ids.begin());

                    if (!pixels_ischeckeds[adjacent_pixel_id]) {
                        active_pixel_ids.push_back(adjacent_pixel_id); // if not checked, then add to the cluster
                        add_adjacent_pixels(adjacent_pixel_id, adjacent_pixel_ids, pixels_ischeckeds, false); // and add its adjacent pixels to the list of adjacent pixels
                        pixels_ischeckeds[adjacent_pixel_id] = true;  // and set it to be checked
                    }
                }

                // calculate energy
                // we should include a layer of calculation to expand the cluster a little bit
                // to include adjacent pixels ADU may not be that large for the clustering, but necessary for the charge sharing
                // since these newly added pixels are directly used for charge sharing correction
                // we put them into the active_pixel_ids
                // for (auto& active_pixel_id : active_pixel_ids)
                //    add_adjacent_pixels(active_pixel_ids, active_pixel_ids, pixels_isvalids, true);
                // pack all active pixels into a cluster
                std::vector<uint16_t> adus;
                std::vector<float> energies;
                float total_energy = 0;

                // // add all neighboring pixels
                if (extended_mode) {
                    std::vector<bool> blank_pixels_ischeckeds(global_config["N_PIXELS"][0], false);
                    for (auto& active_pixel_id : active_pixel_ids) blank_pixels_ischeckeds[active_pixel_id] = true;
                    for (auto& active_pixel_id : active_pixel_ids)
                        add_adjacent_pixels(active_pixel_id, active_pixel_ids, blank_pixels_ischeckeds, true);
                }

                for (auto& active_pixel_id : active_pixel_ids) {
                    adus.push_back(frame_image[active_pixel_id]);
                    float energy = calib0[active_pixel_id] * frame_image[active_pixel_id] + calib1[active_pixel_id];
                    energy = energy > 0 ? energy : 0;
                    total_energy += energy;
                    energies.push_back(energy);
                }
                Cluster cluster = { global_start_frame_id + frame_id, active_pixel_ids, adus, energies};
                threads_clusters[thread_id].push_back(cluster);
            }

            // // print the progress
            // if (frame_id % 5000 == 0 && frame_id != 0) {
            //     auto end_time = std::chrono::high_resolution_clock::now();
            //     auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
            //     int seconds = duration.count();  // This will give you the duration in seconds as an integer
            //     #pragma omp critical
            //     std::cout << format_string("\tThread {} finished frame {} in {} s", thread_id, frame_id, seconds) << std::endl;
            // }
        }

    }
    std::vector<Cluster> clusters;
    for (auto &thread_clusters : threads_clusters)
        clusters.insert(clusters.end(), thread_clusters.begin(), thread_clusters.end());

    return clusters;
}


int save_clusters(const std::string filename, const std::vector<Cluster> &clusters, int pixel_n, bool append) {

    if (!append)  // if not in append mode, then remove the file if it exists
        std::remove(filename.c_str());

    std::ofstream file(filename, std::ios::binary | std::ios::app);
    if (!file.is_open()) throw std::runtime_error("Cannot open: " + filename);

    for (const auto& cluster : clusters) {

        size_t num_pixels = cluster.pixel_ids.size();
        if (pixel_n > 0 && num_pixels != static_cast<size_t>(pixel_n)) continue;
        std::vector<uint16_t> pixel_ids_uint16(cluster.pixel_ids.begin(), cluster.pixel_ids.end()); // this is to make sure the memory alignment, since there are equal number of pixel_ids and adus
        file.write(reinterpret_cast<const char*>(&cluster.frame_id), 4); // int, 4 bytes
        file.write(reinterpret_cast<const char*>(&num_pixels), 4); // int, 4 bytes

        file.write(reinterpret_cast<const char*>(pixel_ids_uint16.data()), num_pixels * sizeof(uint16_t));
        file.write(reinterpret_cast<const char*>(cluster.adus.data()), num_pixels * sizeof(uint16_t));
        file.write(reinterpret_cast<const char*>(cluster.energies.data()), num_pixels * sizeof(float));
    }

    file.close();
    return 0;
}

int raw2clusters(const std::vector<std::string> rawfiles, const std::string crystals_cluster_folder,
    const std::vector<std::string> crystals_calibration_files, const std::vector<std::string> crystals_threshold_files, const bool clear_file, const bool extended_mode) {

    // create output folder if not exists
    if (!std::filesystem::exists(crystals_cluster_folder))
        std::filesystem::create_directory(crystals_cluster_folder);

    std::vector<std::vector<uint8_t>> crystals_pixels_isvalids(global_config["N_CRYSTALS"][0], std::vector<uint8_t>(global_config["N_PIXELS"][0], false));
    std::vector<std::vector<std::vector<float>>> crystals_pixels_peaks(global_config["N_CRYSTALS"][0], std::vector<std::vector<float>>(global_config["N_PIXELS"][0]));
    std::vector<std::vector<std::vector<float>>> crystals_pixels_calibrations(global_config["N_CRYSTALS"][0], std::vector<std::vector<float>>(global_config["N_PIXELS"][0]));
    std::vector<std::vector<float>> crystals_pixels_thresholds(global_config["N_CRYSTALS"][0], std::vector<float>(global_config["N_PIXELS"][0], 0));
    std::vector<std::vector<float>> crystals_pixels_secondary_thresholds(global_config["N_CRYSTALS"][0], std::vector<float>(global_config["N_PIXELS"][0], 0));
    for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
        if (crystals_calibration_files[crystal_id] == "skip") continue;
        std::cout << format_string("Reading calibration and threshold files for crystal {} ...", crystal_id) << std::endl;
        // read pixels calibrations and thresholds
        read_pixels_calibrations(crystals_calibration_files[crystal_id],
            crystals_pixels_calibrations[crystal_id],
            crystals_pixels_peaks[crystal_id],
            crystals_pixels_isvalids[crystal_id],
            crystals_pixels_thresholds[crystal_id]);

        // read pixels thresholds
        if (crystals_threshold_files[crystal_id] != "from_calibration") {
            read_pixels_thresholds(crystals_threshold_files[crystal_id], crystals_pixels_thresholds[crystal_id]);
        }

        // calculate secondary thresholds
        for (int pixel_id = 0; pixel_id < global_config["N_PIXELS"][0]; pixel_id++) {
            crystals_pixels_secondary_thresholds[crystal_id][pixel_id] = global_config["SECONDARY_THRESHOLD"][crystal_id] / crystals_pixels_calibrations[crystal_id][pixel_id][0]; // this is to ensure that thare is at least 1 pixel in the cluster having energy larger than 10 keV + threshold
        }
    }

    // std::cout << "Finished reading calibration and threshold files." << std::endl;
    // std::cout << crystals_pixels_secondary_thresholds[0][0] << std::endl;


    // since the clusters can be very large. we need to use append mode to write them into files.
    for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
        if (crystals_calibration_files[crystal_id] == "skip") continue;
        for (int pixel_n = 0; pixel_n < 10; pixel_n++) {
            std::string cluster_file = format_string("{}\\clusters_crystal_{}_pixel_{}.bin", crystals_cluster_folder, crystal_id, pixel_n);
            if (extended_mode)
                cluster_file = format_string("{}\\clusters_crystal_{}_pixel_{}_extended.bin", crystals_cluster_folder, crystal_id, pixel_n);
            if (clear_file) std::ofstream file(cluster_file, std::ios::binary | std::ios::trunc); // 清空文件
        }
    }

    // define the data container for frames images (flat)
    const int n_crystals = global_config["N_CRYSTALS"][0];
    const int max_frames = global_config["MAX_FRAMES_SIZE"][0];
    const int n_pixels = global_config["N_PIXELS"][0];
    std::vector<uint16_t> crystals_frames_images(static_cast<size_t>(n_crystals) * max_frames * n_pixels, 0);

    // read raw files one by one
    const int n_frame_pixels = global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];
    std::vector<std::vector<Cluster>> crystals_clusters(global_config["N_CRYSTALS"][0]);
    for (auto& rawfile : rawfiles) {

        // read the rawfile by chunks
        std::ifstream file(rawfile, std::ios::binary | std::ios::ate);
        size_t file_size = file.tellg();
        size_t n_frames = file_size / n_frame_pixels / sizeof(uint16_t);
        file.close();

        int n_flow_times = (n_frames + global_config["MAX_FRAMES_SIZE"][0] - 1) / global_config["MAX_FRAMES_SIZE"][0];
        for (int flow_i = 0; flow_i < n_flow_times; flow_i++) {
            int start_frame_id = flow_i * global_config["MAX_FRAMES_SIZE"][0];
            int end_frame_id = start_frame_id + global_config["MAX_FRAMES_SIZE"][0] < n_frames ? start_frame_id + global_config["MAX_FRAMES_SIZE"][0] : n_frames;

            read_rawfile_to_crystals_frames_images(rawfile, start_frame_id, end_frame_id, crystals_frames_images.data());

            for (int crystal_id = 0; crystal_id < n_crystals; crystal_id++) {
                if (crystals_calibration_files[crystal_id] == "skip") continue;
                std::cout << format_string("Start clustering crystal {} from file {}, frames {} to {} ...", crystal_id, rawfile, start_frame_id, end_frame_id) << std::endl;
                const int frames_in_chunk = end_frame_id - start_frame_id;
                const uint16_t* frames_ptr = crystals_frames_images.data() + (static_cast<size_t>(crystal_id) * max_frames) * n_pixels;

                // flatten calibration coefficients for faster access
                std::vector<float> calib0(n_pixels), calib1(n_pixels);
                for (int p = 0; p < n_pixels; ++p) {
                    calib0[p] = crystals_pixels_calibrations[crystal_id][p][0];
                    calib1[p] = crystals_pixels_calibrations[crystal_id][p][1];
                }

                std::vector<Cluster> file_crystal_clusters = read_frames_images_to_clusters(
                    frames_ptr, frames_in_chunk, calib0.data(), calib1.data(),
                    crystals_pixels_isvalids[crystal_id].data(), crystals_pixels_thresholds[crystal_id].data(),
                    crystals_pixels_secondary_thresholds[crystal_id].data(), start_frame_id, extended_mode);
                crystals_clusters[crystal_id].insert(crystals_clusters[crystal_id].end(), file_crystal_clusters.begin(), file_crystal_clusters.end());
            }

            for (int crystal_id = 0; crystal_id < global_config["N_CRYSTALS"][0]; crystal_id++) {
                if (crystals_calibration_files[crystal_id] == "skip") continue;
                std::string cluster_file = format_string("{}\\clusters_crystal_{}.bin", crystals_cluster_folder, crystal_id);
                for (int pixel_n = 1; pixel_n <= 9; pixel_n++) {
                    std::string cluster_file = format_string("{}\\clusters_crystal_{}_pixel_{}.bin", crystals_cluster_folder, crystal_id, pixel_n);
                    if (extended_mode)
                        cluster_file = format_string("{}\\clusters_crystal_{}_pixel_{}_extended.bin", crystals_cluster_folder, crystal_id, pixel_n);
                    save_clusters(cluster_file, crystals_clusters[crystal_id], pixel_n, true);
                }
                crystals_clusters[crystal_id].clear();
            }

        }
    }
    return 0;
}
