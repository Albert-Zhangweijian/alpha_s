#include <vector>
#include <string>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <numeric>
#include <windows.h>

#include "header.hpp"



uint16_t* mmap_windows(const std::string& path, size_t& file_size, HANDLE& hMap)
{
    HANDLE hFile = CreateFileA(path.c_str(),
        GENERIC_READ, FILE_SHARE_READ,
        NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    if (hFile == INVALID_HANDLE_VALUE) {
        return nullptr;
    }

    LARGE_INTEGER size;
    if (!GetFileSizeEx(hFile, &size)) {
        CloseHandle(hFile);
        return nullptr;
    }

    file_size = static_cast<size_t>(size.QuadPart);   // 64-bit OK

    hMap = CreateFileMappingA(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    CloseHandle(hFile);

    if (!hMap) {
        return nullptr;
    }

    return (uint16_t*)MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, 0);
}

int read_single_pixel_batch(const std::vector<std::string>& rawfiles, const std::vector<int>& crystal_ids, const std::vector<int>& rows, const std::vector<int>& cols, std::vector<std::vector<uint16_t>>& multiple_pixel_values) {

    static std::vector<int> pixel_order_vec(global_config["N_PIXELS_PREMERGE"][0]);
    compute_pixel_order(
        global_config["N_COLS_PREMERGE"][0], global_config["N_ROWS_PREMERGE"][0],
        global_config["N_READOUT_PIXELS"][0], global_config["N_READOUT_GROUPS"][0],
        global_config["N_PIXELS_PREMERGE"][0], global_config["N_PIXELS"][0],
        global_config["ROW_MERGE_FOLD"][0], global_config["COL_MERGE_FOLD"][0],
        global_config["ROW_MERGE_INDEX"][0], global_config["COL_MERGE_INDEX"][0],
        global_config["N_COLS"][0],
        pixel_order_vec.data()
    );

    // only the position of the pixel in the pixel array
    int n_pixels = crystal_ids.size();
    std::vector<int> pixel_position_offsets(n_pixels, 0);  // the offset of this pixel in the flattened pixel array
    for (int i = 0; i < n_pixels; ++i) {
        int pixel_id = rows[i] * global_config["N_COLS"][0] + cols[i];
        for (int index = 0; index < global_config["N_PIXELS_PREMERGE"][0]; ++index) {
            if (pixel_order_vec[index] == pixel_id) {
                pixel_position_offsets[i] = index;
                break;
            }
        }
    }

    // cause every crystal's same position pixel is stored together
    // the offset in every single file
    for (int i = 0; i < n_pixels; ++i) {
        pixel_position_offsets[i] = pixel_position_offsets[i] * global_config["N_CRYSTALS"][0] + crystal_ids[i]; 
    }

    // calculate total number of frames
    std::vector<size_t> n_frames_per_file;
    int offset_per_frame = global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];
    for (const auto& rawfile : rawfiles) {
        std::ifstream file(rawfile, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << rawfile << std::endl;
            continue;
        }
        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        size_t n_frames_in_file = file_size / offset_per_frame / sizeof(uint16_t);
        n_frames_per_file.push_back(n_frames_in_file);
        file.close();
    }
    size_t n_frames_total = std::accumulate(n_frames_per_file.begin(), n_frames_per_file.end(), 0ULL);

    // resize the data structure to hold all pixel values
    multiple_pixel_values.assign(n_pixels, std::vector<uint16_t>(n_frames_total));

    // read each pixel one by one
    int file_frame_start_index = 0;
    uint16_t pixel_value;
    for (int file_i = 0; file_i < rawfiles.size(); ++file_i) {
        const auto& rawfile = rawfiles[file_i];
        std::cout << "Reading file " << file_i + 1 << " / " << rawfiles.size() << ": " << rawfile << "\n";

        size_t file_size_get;
        HANDLE hMap;
        uint16_t* file_ptr = mmap_windows(rawfile, file_size_get, hMap);
        if (file_size_get != n_frames_per_file[file_i] * offset_per_frame * sizeof(uint16_t)) {
            std::cout << "sizeof(void*) = " << sizeof(void*) << "\n";
            std::cerr << "Error: file size mismatch for file: " << rawfile << std::endl;
            std::cerr << "Expected size: " << n_frames_per_file[file_i] * offset_per_frame << ", got: " << file_size_get << std::endl;
            UnmapViewOfFile(file_ptr);
            CloseHandle(hMap);
            continue;
        }
         
        for (size_t frame_id = 0; frame_id < n_frames_per_file[file_i]; ++frame_id) {
            file_ptr += offset_per_frame;  // move to the next frame start
            for (int pixel_i = 0; pixel_i < n_pixels; ++pixel_i) {
                multiple_pixel_values[pixel_i][file_frame_start_index + frame_id] = *(file_ptr + pixel_position_offsets[pixel_i]);
            }
        }

        UnmapViewOfFile(file_ptr);
        CloseHandle(hMap);
        file_frame_start_index += n_frames_per_file[file_i];
        std::cout << "Finished reading file " << file_i + 1 << " / " << rawfiles.size() << "\n";
    }

    return 0;
}

int instant_read_batch(const std::vector<std::string>& rawfiles, const std::vector<int>& crystal_ids, const std::vector<int>& rows, const std::vector<int>& cols, std::string type, std::string output_folder) {

    std::vector<std::vector<uint16_t>> pixel_values;
    read_single_pixel_batch(rawfiles, crystal_ids, rows, cols, pixel_values);

    if (type == "scatter") {
        // compute and save scatter plot
        for (size_t i = 0; i < crystal_ids.size(); ++i) {
            int crystal_id = crystal_ids[i];
            int row = rows[i];
            int col = cols[i];
            std::string output_path = output_folder + "/crystal_" + std::to_string(crystal_id) + "_pixel_" + std::to_string(row) + "_" + std::to_string(col) + "_scatter.bin";
            std::ofstream outfile(output_path, std::ios::binary);
            outfile.write(reinterpret_cast<const char*>(pixel_values[i].data()), pixel_values[i].size() * sizeof(uint16_t));
            outfile.close();
            std::cout << "Scatter data saved to " << output_path << "\n";
        }
    } else if (type == "spectrum") {
        // compute and save spectrum plot
        std::vector<float> bins = create_bins(global_config["N_BINS"][0], global_config["ADU_MIN"][0], global_config["ADU_MAX"][0]);
        static const std::vector<size_t> lookup_table = compute_adu_lookup_table(bins);

        std::vector<uint32_t> spectrum(global_config["N_BINS"][0], 0);
        for (const auto& values : pixel_values) {
            for (const auto& value : values) {
                spectrum[lookup_table[static_cast<size_t>(value - global_config["ADU_MIN"][0])]]++;
            }
        }

        for (size_t i = 0; i < crystal_ids.size(); ++i) {
            int crystal_id = crystal_ids[i];
            int row = rows[i];
            int col = cols[i];
            std::string output_path = output_folder + "/crystal_" + std::to_string(crystal_id) + "_pixel_" + std::to_string(row) + "_" + std::to_string(col) + "_spectrum.bin";
            std::ofstream outfile(output_path, std::ios::binary);
            outfile.write(reinterpret_cast<const char*>(spectrum.data()), spectrum.size() * sizeof(uint32_t));
            outfile.close();
            std::cout << "Spectrum data saved to " << output_path << "\n";
        }
    }
    return 0;
}