#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <vector>

#include "header.hpp"

std::vector<int> split(const std::string& line) {
    std::istringstream iss(line);
    std::vector<int> tokens;
    int token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::unordered_map<std::string, std::vector<int>>  read_config(const std::string& filename, bool verbose) {
    std::unordered_map<std::string, std::vector<int>> config;
    if (!std::filesystem::exists(filename))
        throw std::invalid_argument("Config file does not exist: " + filename);
    if (std::filesystem::path(filename).extension() != ".txt")
        throw std::invalid_argument("Config file must be a .txt file");
    std::ifstream infile(filename);
    if (!infile.is_open()) throw std::runtime_error("Could not open config file: " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        // neglect comments
        auto pos_comment = line.find_first_of("#/");
        if (pos_comment != std::string::npos) {
            line = line.substr(0, pos_comment);
        }

        // remove leading and trailing whitespace
        if (line.empty()) continue;

        // find key=value
        auto pos = line.find('=');
        if (pos == std::string::npos) continue; // invalid line

        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);

        // 去掉多余空格
        auto trim = [](std::string& s) {
            size_t start = s.find_first_not_of(" \t");
            size_t end = s.find_last_not_of(" \t");
            if (start == std::string::npos) { s = ""; return; }
            s = s.substr(start, end - start + 1);
        };
        trim(key);
        trim(value);

        // 拆分 value（支持多个值）
        std::vector<int> values = split(value);
        config[key] = values;
    }

    if (verbose) {
        std::cout << "Configuration loaded from " << filename << ":\n";
        for (const auto& [key, values] : config) {
            std::cout << key << " = ";
            for (const auto& val : values) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    return config;
}


std::unordered_map<std::string, std::vector<int>> global_config;