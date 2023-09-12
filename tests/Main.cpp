#include <catch2/catch_session.hpp>
#include <iostream>

#include "TestData.h"

int main(int argc, char* argv[]) {
    Catch::Session session; // There must be exactly one instance

    std::string data_path; // Some user variable you want to be able to set

    // Build a new parser on top of Catch2's
    using namespace Catch::Clara;
    auto cli
        = session.cli()           // Get Catch2's command line parser
        | Opt(data_path, "data path")["--data_path"]("test cases data path");        // description string for the help output

    // Now pass the new composite back to Catch2 so it uses that
    session.cli(cli);

    // Let Catch2 (using Clara) parse the command line
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

    // if set on the command line then 'height' is now set at this point
    if (!data_path.empty()) {
        get_test_data(data_path);
    }

    if (TEST_DATA_PP.empty() && TEST_DATA_PM.empty() && TEST_DATA_MP.empty() && TEST_DATA_MM.empty()) {
        std::cout << "No test data found. Please set the data path with -d or --data_path" << std::endl;
        return -1;
    }

    return session.run();
}
