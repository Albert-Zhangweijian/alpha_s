#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <cstdint>
#include <utility>
#include <unordered_map>

// config
extern std::unordered_map<std::string, std::vector<int>> global_config;
std::unordered_map<std::string, std::vector<int>> read_config(const std::string& filename, bool verbose);


// basics
std::vector<float> create_bins(unsigned int n_bins, float low, float high);
int read_pixels_thresholds(const std::string filename, std::vector<float>& pixels_thresholds);
int read_pixels_calibrations(const std::string filename,
    std::vector<std::vector<float>>& pixels_calibrations,
    std::vector<std::vector<float>>& pixels_peaks,
    std::vector<uint8_t>& pixels_isvalids,
    std::vector<float>& pixels_thresholds);
int read_rawfile_to_crystals_frames_images(const std::string& rawfile, const int start_frame_id, const int end_frame_id, uint16_t* crystals_frames_images);

// functional modules
int raw2frames(const std::vector<std::string> rawfiles, const std::string crystals_frames_folder);
int raw2spectra(const std::vector<std::string> rawfiles, const std::string crystals_spectra_folder);
int raw2scatters(const std::vector<std::string>& rawfiles, const std::string crystals_scatters_folder,
    const std::vector<std::string>& crystals_calibration_files, const std::string process_type);
int raw2clusters(const std::vector<std::string> rawfiles, const std::string crystals_cluster_folder,
    const std::vector<std::string> crystals_calibration_files, const std::vector<std::string> crystals_threshold_files,
    const bool plot_frames_and_prevent_clear, const bool extended_mode);

// string format
// template <typename... Args>
// std::string format_string(const std::string& string_format, Args&&... args) {
//     std::vector<std::string> arg_list{std::to_string(std::forward<Args>(args))...};
//     std::ostringstream ostream;
//     size_t arg_index = 0;

//     for (size_t i = 0; i < string_format.size(); ++i) {
//         if (string_format[i] == '{' && i + 1 < string_format.size() && string_format[i + 1] == '}') {
//             if (arg_index < arg_list.size()) {
//                 ostream << arg_list[arg_index++];
//             } else {
//                 ostream << "{}";  // if there are not enough arguments, keep the placeholder
//             }
//             ++i;  // skip the next '}'
//         } else {
//             ostream << string_format[i];
//         }
//     }
//     return ostream.str();
// }
template <typename... Args>
std::string format_string(const std::string& string_format, Args&&... args) {
    // 内部lambda，把任意类型转成 std::string
    auto to_string_any = [](auto&& value) -> std::string {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, std::string>) {
            return value;
        } else if constexpr (std::is_same_v<T, const char*>) {
            return std::string(value);
        } else {
            std::ostringstream oss;
            oss << value;
            return oss.str();
        }
    };

    std::vector<std::string> arg_list{to_string_any(std::forward<Args>(args))...};
    std::ostringstream ostream;
    size_t arg_index = 0;

    for (size_t i = 0; i < string_format.size(); ++i) {
        if (string_format[i] == '{' && i + 1 < string_format.size() && string_format[i + 1] == '}') {
            if (arg_index < arg_list.size()) {
                ostream << arg_list[arg_index++];
            } else {
                ostream << "{}";  // 没有足够参数就保留
            }
            ++i;  // 跳过 '}'
        } else {
            ostream << string_format[i];
        }
    }
    return ostream.str();
}
