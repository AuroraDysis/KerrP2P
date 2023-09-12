#include "TestData.h"

#include <string>

#include <boost/filesystem.hpp>

void read_csv(std::vector<std::array<std::string, 9>>& data, const boost::filesystem::path& file_path) {
	std::ifstream ifs(file_path.string());
	if (!ifs.is_open()) {
		return;
	}
	std::string line;
  std::array<std::string, 9> row;
	while (std::getline(ifs, line)) {
		std::stringstream ss(line);
		std::string cell;
		int i = 0;
		while (std::getline(ss, cell, ',')) {
			row[i++] = cell;
		}
    data.emplace_back(row);
	}
}

void get_test_data(std::string& path) {
    using namespace boost::filesystem;

    if (!is_directory(path) || !exists(path)) {
        return;
    }

    read_csv(TEST_DATA_PP, path / "pp.csv");
    read_csv(TEST_DATA_PM, path / "pm.csv");
    read_csv(TEST_DATA_MP, path / "mp.csv");
    read_csv(TEST_DATA_MM, path / "mm.csv");
}
