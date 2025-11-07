#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "cxxopts.hpp"
#include "header.hpp"

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <command> [options]\n";
        std::cout << "Available commands: raw2spectra, raw2scatters, and raw2clusters\n";
        return 0;
    }

    std::string command = argv[1];  // the first is the subcommand
    std::vector<char*> args(argv + 2, argv + argc); // the remaining arguments
    std::cout << "Executing command: " << command << "\n" << std::endl;
    
    global_config = read_config("config.txt", false);
    


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

        options.add_options()
            ("i,input", "Raw input files", cxxopts::value<std::vector<std::string>>(rawfiles))
            ("c,calibration", "Use calibration files", cxxopts::value<std::vector<std::string>>(crystals_calibration_files))
            ("t,threshold", "Use threshold files", cxxopts::value<std::vector<std::string>>(crystals_threshold_files))
            ("o,output", "Output folder", cxxopts::value<std::string>(crystals_cluster_folder)->default_value("clusters"))
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


        raw2clusters(rawfiles, crystals_cluster_folder,
            crystals_calibration_files, crystals_threshold_files, true);

        return 0;
    }

    std::cout << "Unknown command: " << command << "\n";
    std::cout << "Available commands: raw2spectra, raw2scatters, and raw2clusters\n";
    return 1;

}
