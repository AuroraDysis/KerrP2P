#include "TestData.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace boost::property_tree;

void get_test_data(std::string &path) {
  using namespace boost::filesystem;
  std::vector<ptree> test_data;

  if (!is_directory(path) || !exists(path)) {
    return;
  }

  for(auto& entry : boost::make_iterator_range(directory_iterator(path), {})) {
    // find .json files and parse them
    if (entry.path().extension() == ".json") {
      std::ifstream file(entry.path().string());
      ptree jsontree;
      read_json(file, jsontree);
      TEST_DATA.push_back(std::move(jsontree));
    }
  }
}
