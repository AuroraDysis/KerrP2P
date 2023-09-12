#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "ForwardRayTracing.h"
#include "Utils.h"

using std::string;

int main(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("a", po::value<std::string>(), "black hole spin");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("compression")) {
        std::cout << "Compression level was set to "
                  << vm["compression"].as<int>() << ".\n";
    } else {
        std::cout << "Compression level was not set.\n";
    }
}
