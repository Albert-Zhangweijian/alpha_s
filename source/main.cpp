#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "cxxopts.hpp"
#include "header.hpp"

int main(int argc, char* argv[]) {

    if (argc <= 2) {
        std::cout << "Usage: " << argv[0] << " <command> [options]\n";
        std::cout << "Available commands: raw2spectra, raw2scatters, and raw2clusters\n";
        return 0;
    }

    // Check for config argument
    std::string config_filepath = "config.txt";
    if (std::string(argv[1]) == "--config")
        config_filepath = argv[2];
    global_config = read_config(config_filepath, false);
    
    // Shift arguments if config is provided
    int shift = (std::string(argv[1]) == "--config") ? 2 : 0;
    std::string command = argv[1 + shift];  // the first is the subcommand
    std::vector<char*> args(argv + 2 + shift, argv + argc); // the remaining arguments
    std::cout << "Executing command: " << command << " with config file: " << config_filepath << "\n" << std::endl;

    if (command == "raw2frames") {

        cxxopts::Options options("alpha raw2frames", "Convert rawfiles to frames");
        std::vector<std::string> rawfiles;
        std::string crystals_frames_folder;

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("o,output", "Output folder", cxxopts::value<std::string>(crystals_frames_folder)->default_value("frames"))
            ("h,help", "Print usage");
        options.parse_positional({"input"});
        auto result = options.parse(args.size(), args.data());
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
        for (const auto& file : rawfiles)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Output: " << crystals_frames_folder << "\n";

        raw2frames(rawfiles, crystals_frames_folder);

        return 0;
    }

    if (command == "raw2spectra") {

        cxxopts::Options options("alpha raw2spectra", "Convert rawfiles to spectra");
        std::vector<std::string> rawfiles;
        std::string crystals_spectra_folder;

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("o,output", "Output folder", cxxopts::value<std::string>(crystals_spectra_folder)->default_value("spectra"))
            ("h,help", "Print usage");
        options.parse_positional({"input"});
        auto result = options.parse(args.size(), args.data());
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
        for (const auto& file : rawfiles)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Output: " << crystals_spectra_folder << "\n";

        raw2spectra(rawfiles, crystals_spectra_folder);

        return 0;
    } 
    
    
    if (command == "raw2scatters") {

        cxxopts::Options options("alpha raw2scatters", "Convert rawfiles to scatters");
        std::vector<std::string> rawfiles;
        std::vector<std::string> crystals_calibration_files;
        std::string crystals_scatters_folder;
        std::string mode;

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("c,calibration", "Use calibration files", cxxopts::value<std::vector<std::string>>(crystals_calibration_files))
            ("o,output", "Output folder", cxxopts::value<std::string>(crystals_scatters_folder)->default_value("scatters"))
            ("m,mode", "Process mode", cxxopts::value<std::string>(mode)->default_value("random"))
            ("h,help", "Print usage");
        options.parse_positional({"input", "calibration"});
        auto result = options.parse(args.size(), args.data());
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }


        std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
        for (const auto& file : rawfiles)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Found " << crystals_calibration_files.size() << " calibration files:" << std::endl;
        for (const auto& file : crystals_calibration_files)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Output: " << crystals_scatters_folder << "\n";
        std::cout << "Mode: " << mode << "\n";


        raw2scatters(rawfiles, crystals_scatters_folder, crystals_calibration_files, mode);
        return 0;
    }

    if (command == "raw2clusters") {

        cxxopts::Options options("alpha raw2clusters", "Convert rawfiles to clusters");
        std::vector<std::string> rawfiles;
        std::vector<std::string> crystals_calibration_files;
        std::vector<std::string> crystals_threshold_files;
        std::string crystals_cluster_folder;
        bool extended_mode = false;

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("c,calibration", "Use calibration files", cxxopts::value<std::vector<std::string>>(crystals_calibration_files))
            ("t,threshold", "Use threshold files", cxxopts::value<std::vector<std::string>>(crystals_threshold_files))
            ("o,output", "Output folder", cxxopts::value<std::string>(crystals_cluster_folder)->default_value("clusters"))
            ("e,extended_mode", "Enable output exteneded clusters", cxxopts::value<bool>(extended_mode)->default_value("false")->implicit_value("true"))
            ("h,help", "Print usage");
        options.parse_positional({"input", "calibration", "threshold"});
        auto result = options.parse(args.size(), args.data());
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
        for (const auto& file : rawfiles)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Found " << crystals_calibration_files.size() << " calibration files:" << std::endl;
        for (const auto& file : crystals_calibration_files)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Found " << crystals_threshold_files.size() << " threshold files:" << std::endl;
        for (const auto& file : crystals_threshold_files)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Output: " << crystals_cluster_folder << "\n";
        std::cout << "Extended mode: " << (extended_mode ? "true" : "false") << "\n";


        raw2clusters(rawfiles, crystals_cluster_folder,
            crystals_calibration_files, crystals_threshold_files, true, extended_mode);

        return 0;
    }

    // if (command == "instant") {

    //     cxxopts::Options options("alpha instant", "Instant read from rawfiles");
    //     std::vector<std::string> rawfiles;
    //     int crystal_id;
    //     int row;
    //     int col;
    //     std::string type;
    //     std::string output_folder;

    //     options.add_options()
    //         ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
    //         ("crystal_id", "Crystal ID", cxxopts::value<int>(crystal_id))
    //         ("row", "Row index", cxxopts::value<int>(row))
    //         ("col", "Column index", cxxopts::value<int>(col))
    //         ("t,type", "Data type (e.g., 'spectrum', 'scatter')", cxxopts::value<std::string>(type))
    //         ("o,output", "Output folder", cxxopts::value<std::string>(output_folder)->default_value("instant_output"))
    //         ("h,help", "Print usage");
    //     options.parse_positional({"input", "crystal_id", "row", "col"});
    //     auto result = options.parse(args.size(), args.data());
    //     if (result.count("help")) {
    //         std::cout << options.help() << std::endl;
    //         return 0;
    //     }

    //     std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
    //     for (const auto& file : rawfiles)
    //         std::cout << "- [" << file << "]" << std::endl;
    //     std::cout << "Crystal ID: " << crystal_id << "\n";
    //     std::cout << "Row: " << row << "\n";
    //     std::cout << "Column: " << col << "\n";
    //     std::cout << "Type: " << type << "\n";
    //     std::cout << "Output folder: " << output_folder << "\n";

    //     instant_read(rawfiles, crystal_id, row, col, type, output_folder);

    //     return 0;
    // } 

    if (command == "instant") {

        cxxopts::Options options("alpha instant", "Instant read from rawfiles");
        std::vector<std::string> rawfiles;
        std::vector<int> crystal_ids;
        std::vector<int> rows;
        std::vector<int> cols;
        std::string type;
        std::string output_folder;

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("crystal_ids", "Crystal IDs", cxxopts::value<std::vector<int>>(crystal_ids))
            ("rows", "Row indices", cxxopts::value<std::vector<int>>(rows))
            ("cols", "Column indices", cxxopts::value<std::vector<int>>(cols))
            ("t,type", "Data type (e.g., 'spectrum', 'scatter')", cxxopts::value<std::string>(type))
            ("o,output", "Output folder", cxxopts::value<std::string>(output_folder)->default_value("instant_output"))
            ("h,help", "Print usage");
        options.parse_positional({"input", "crystal_ids", "rows", "cols"});
        auto result = options.parse(args.size(), args.data());
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        std::cout << "Found " << rawfiles.size() << " raw files:" << std::endl;
        for (const auto& file : rawfiles)
            std::cout << "- [" << file << "]" << std::endl;
        std::cout << "Crystal IDs: ";
        for (const auto& id : crystal_ids) std::cout << id << " ";
        std::cout << "\n";
        std::cout << "Rows: ";
        for (const auto& r : rows) std::cout << r << " ";
        std::cout << "\n";
        std::cout << "Columns: ";
        for (const auto& c : cols) std::cout << c << " ";
        std::cout << "\n";
        std::cout << "Type: " << type << "\n";
        std::cout << "Output folder: " << output_folder << "\n";

        instant_read_batch(rawfiles, crystal_ids, rows, cols, type, output_folder);

        return 0;
    }    

    std::cout << "Unknown command: " << command << "\n";
    std::cout << "Available commands: raw2spectra, raw2scatters, and raw2clusters\n";
    return 1;

}
