#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <string.h>
#include "header.hpp"

// calculate the index order
// the raw data is in order like this:
// (Frame 1) (Frame 2) (Frame 3) ... (Frame 10000)
// Frame 1 = (Pixel 1,1) (Pixel 1,21) (Pixel 1,41) (Pixel 1,61) (Pixel 2,1) (Pixel 2,21) ...
// Pixel 1,1 = (Crystal 1-1,1) (Crystal 2-1,1) (Crystal 3-1,1) (Crystal 4-1,1)
// col_order is the real column of each pixel in one row
// pixel_order is the order of each pixel in one frame
int compute_pixel_order(int n_cols_premerge, int n_rows_premerge, int n_readout_pixels, int n_readout_groups, int n_pixels_premerge, int n_pixels,
    int row_merge_fold, int col_merge_fold, int row_merge_index, int col_merge_index, int n_cols, int* pixel_order) {

    // compute the column order in a row
    std::vector<int> col_order(n_cols_premerge);
    int index = 0;
    for (int pin_id = 0; pin_id < n_readout_pixels; ++pin_id) {
        for (int band_id = 0; band_id < n_readout_groups; ++band_id) {
            col_order[index] = pin_id + band_id * n_readout_pixels;
            index++;
        }
    }

    // compute the pixel order in a crystal
    //the data inside raw files will follow this order
    index = 0;
    for (int row = 0; row < n_rows_premerge; ++row) {
        for (int col_i = 0; col_i < n_cols_premerge; ++col_i) {
            int col = col_order[col_i];
            pixel_order[index] = row * n_cols_premerge + col;
            index++;
        }
    }


    // considering the pixel merging technique
    // map index pre-merge to index post-merge
    // the vector size is still the same, but pixels without data will be marked as -1
    if (n_pixels != n_pixels_premerge) {

        std::vector<int> pixel_order_postmerge(n_pixels_premerge);
        for (int index = 0; index < n_pixels_premerge; ++index) {
            int pixel_id_premerge = pixel_order[index];
            int row_premerge = pixel_id_premerge / n_cols_premerge;
            int col_premerge = pixel_id_premerge % n_cols_premerge;

            if (row_premerge % row_merge_fold != row_merge_index || col_premerge % col_merge_fold != col_merge_index) {
                pixel_order_postmerge[index] =  -1;
                continue; // skip the pixels that are not the first pixel in the merged group
            }
            int row_postmerge = (row_premerge - row_merge_index) / row_merge_fold;
            int col_postmerge = (col_premerge - col_merge_index) / col_merge_fold;
            int pixel_id_postmerge = row_postmerge * n_cols + col_postmerge;
            pixel_order_postmerge[index] = pixel_id_postmerge;
        }

        // update the pixel order to the merged one
        for (int index = 0; index < n_pixels_premerge; ++index)
            pixel_order[index] = pixel_order_postmerge[index];
    }

    return 0;
}

int read_rawfile_to_crystals_frames_images(const std::string& rawfile, const int start_frame_id, const int end_frame_id, uint16_t* crystals_frames_images) {

    auto start_time = std::chrono::high_resolution_clock::now();

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

    // retrive the only filename from the path
    size_t last_slash_idx = rawfile.find_last_of("\\/");
    std::string filename;
    if (std::string::npos != last_slash_idx) {
        filename = rawfile.substr(last_slash_idx + 1);
    } else {
        filename = rawfile;
    }

    // 预分配三维存储结构
    int n_frames_to_read = end_frame_id - start_frame_id;
    const size_t n_frame_pixels = global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];

    // 分块处理参数设置
    const size_t n_chunk_frames = 2000;
    std::vector<uint16_t> read_buffer(n_chunk_frames * n_frame_pixels);

    // 主处理循环
    std::ifstream file(rawfile, std::ios::binary);
    file.seekg(start_frame_id * n_frame_pixels * sizeof(uint16_t));
    for (size_t n_readed_frames = 0; n_readed_frames < n_frames_to_read; n_readed_frames += n_chunk_frames) {
        const size_t n_current_chunk_frames = std::min(n_chunk_frames, n_frames_to_read - n_readed_frames);
        const size_t bytes_to_read = n_current_chunk_frames * n_frame_pixels * sizeof(uint16_t);

        // 读取当前块数据
        file.read(reinterpret_cast<char*>(read_buffer.data()), bytes_to_read);

        // 并行处理当前块
        #pragma omp parallel for num_threads(global_config["N_THREADS"][0])
        for (int current_frame_i = 0; current_frame_i < n_current_chunk_frames; ++current_frame_i) {
            const uint16_t* frame_start = read_buffer.data() + current_frame_i * global_config["N_CRYSTALS"][0] * global_config["N_PIXELS_PREMERGE"][0];
            const size_t global_frame_i = n_readed_frames + current_frame_i;

            // 优化后的像素处理循环
            for (int crystal_i = 0; crystal_i < global_config["N_CRYSTALS"][0]; ++crystal_i) {
                const size_t frame_offset = (static_cast<size_t>(crystal_i) * global_config["MAX_FRAMES_SIZE"][0] + global_frame_i) * static_cast<size_t>(global_config["N_PIXELS"][0]);
                uint16_t* target_frame = crystals_frames_images + frame_offset;

                // 使用预计算索引加速访问
                for (int pixel_preorder_i = 0; pixel_preorder_i < global_config["N_PIXELS_PREMERGE"][0]; ++pixel_preorder_i) {
                    int mapped = pixel_order_vec[pixel_preorder_i];
                    if (mapped != -1) {
                        target_frame[mapped] = frame_start[pixel_preorder_i * global_config["N_CRYSTALS"][0] + crystal_i];
                        // std::cout << format_string("\rCrystal {}, Frame {}, Pixel {}: Value {}", crystal_i, global_frame_i, pixel_order[pixel_preorder_i], target_frame[pixel_order[pixel_preorder_i]]) << std::endl;
                    }
                }
            }
        }

        // 进度显示
        std::cout << format_string("\rReading {} frames ({}%) from file {}", n_readed_frames + n_current_chunk_frames, (n_readed_frames + n_current_chunk_frames)*100.0/n_frames_to_read, filename) << std::flush;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << format_string("\nFinished {} frames ({}%) from file {},  {} s elapsed. IO speed = {} Mb/s", n_frames_to_read, (n_frames_to_read)*100.0/n_frames_to_read, filename, elapsed.count(), (n_frames_to_read * n_frame_pixels * sizeof(uint16_t)) / elapsed.count() / (1024 * 1024)) << std::endl;
    file.close();

    // set crystals_frames_images to zero if not all frames are read
    if (end_frame_id - start_frame_id != n_frames_to_read) {
        size_t start_zero = static_cast<size_t>(end_frame_id - start_frame_id) * global_config["N_PIXELS"][0];
        size_t zero_count = static_cast<size_t>(n_frames_to_read - (end_frame_id - start_frame_id)) * global_config["N_PIXELS"][0];
        memset(crystals_frames_images + start_zero, 0, zero_count * sizeof(uint16_t));
    }


    return 0;
}

