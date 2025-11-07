#include <vector>
#include <string>
#include <stdexcept>
#include "highfive/HighFive.hpp"

#include "header.hpp"

std::vector<float> create_bins(unsigned int n_bins, float low, float high) {
    if (low >= high) throw std::invalid_argument("E_min must be less than E_max");
    float step = (high - low) / n_bins;
    std::vector<float> bins(n_bins);
    for (unsigned int i = 0; i < n_bins; i++) bins[i] = low + i * step;
    return bins;
}

int read_pixels_thresholds(const std::string filename, std::vector<float>& pixels_thresholds) {
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getDataSet("pixels_thresholds").read(pixels_thresholds);
    return 0;
}

int read_pixels_calibrations(const std::string filename,
    std::vector<std::vector<float>>& pixels_calibrations,
    std::vector<std::vector<float>>& pixels_peaks,
    std::vector<uint8_t>& pixels_isvalids,
    std::vector<float>& pixels_thresholds) {
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getDataSet("pixels_calibrations").read(pixels_calibrations);
    // optional datasets
    if (file.exist("pixels_peaks")) file.getDataSet("pixels_peaks").read(pixels_peaks);
    if (file.exist("pixels_isvalids")) file.getDataSet("pixels_isvalids").read(pixels_isvalids);
    if (file.exist("pixels_thresholds")) file.getDataSet("pixels_thresholds").read(pixels_thresholds);
    return 0;
}
