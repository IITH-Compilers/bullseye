// /*
//  * Copyright (c) 2022, BullsEye
//    Adopted from HayStack (PLDI 2019)
// */

// #include<iostream>
// #include"lib/bullseyelib.h"
// int main() {
//   std::cout <<" Hey this is the main function" << std::endl;
//   bullseye();
//   return 0;
// }

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

#include <isl/options.h>

#include "lib/bullseyelib.h"

namespace po = boost::program_options;

const int CACHE_LINE_SIZE = 64;

// default
const int CACHE_SIZE1 = 32 * 1024;
const int CACHE_SIZE2 = 512 * 1024;

// define print operators
namespace std {
std::ostream &operator<<(std::ostream &os, const std::vector<long> &vec) {
  for (int i = 0; i < vec.size(); ++i) {
    os << vec[i];
    if (i < vec.size() - 1)
      os << " ";
  }
  return os;
}
} // namespace std

bullseyelib::ProgramParameters
ConvertToProgramParameters(po::variables_map Variables) {
  std::string InputFile = Variables["input-file"].as<std::string>();
  bullseyelib::ProgramParameters PP(InputFile);
  PP.CacheSizes =
      Variables["cache-sizes"].as<std::vector<long>>();
  PP.LineSize = Variables["line-size"].as<long>();
  
  if (Variables.count("include-path") > 0) {
    std::vector<std::string> IncludePath =
        Variables["include-path"].as<std::vector<std::string>>();
    PP.IncludePath = IncludePath;
  }
  if (Variables.count("define-parameters") > 0) {
    PP.DefineParameters =
        Variables["define-parameters"].as<std::vector<std::string>>();
  }
  if (Variables.count("scop-function")) {
    PP.ScopFunction = Variables["scop-function"].as<std::string>();
  }
  if (Variables.count("compute-counds")) {
    PP.ComputeBounds = Variables["compute-bounds"].as<bool>();
  }
  return PP;
}


int main(int argc, const char **args) {
  try {
    // define the program options
    po::options_description Descriptor("Program options");
    Descriptor.add_options()                    //
        ("help,h", "print the program options") //
        ("cache-sizes,c",
         po::value<std::vector<long>>()->multitoken()->default_value(
             {CACHE_SIZE1, CACHE_SIZE2}),
         "cache sizes in byte") //
        ("line-size,l", po::value<long>()->default_value(CACHE_LINE_SIZE),
         "cache-line size in byte") //
        ("input-file,f", po::value<std::string>(),
         "set the source file [file name]") //
        ("include-path,I", po::value<std::vector<std::string>>(),
         "set the include path [include path]") //
        ("define-parameters,d",
         po::value<std::vector<std::string>>()->multitoken(),
         "parameter values [N=10 M=100]") //
        ("scop-function,s", po::value<std::string>(),
         "set the scop function scop") //
        ("compute-bounds,b", po::value<bool>()->default_value(false),
         "compute stack distance bounds");

    // parse the program options
    po::variables_map Variables;
    po::store(po::parse_command_line(argc, args, Descriptor), Variables);
    po::notify(Variables);
    if (Variables.count("help") || Variables.count("input-file") == 0) {
      std::cout << Descriptor << std::endl;
      return 0;
    }
    // check if the include paths are valid
    for (int i = 0; i < Variables.count("include-path"); ++i) {
      std::string IncludePath =
          Variables["include-path"].as<std::vector<std::string>>()[i];
      if (!bullseyelib::check_path(IncludePath)) {
        printf("-> exit(-1) include path %s not valid\n", IncludePath.c_str());
        exit(-1);
      }
    }
    // check if the source file is valid
    if (!bullseyelib::check_path(Variables["input-file"].as<std::string>())) {
      printf("-> exit(-1) input file %s not found\n",
             Variables["input-file"].as<std::string>().c_str());
      exit(-1);
    }
    // allocate the context outside of the cache model
    std::vector<std::string> IncludePaths;
    if (Variables.count("include-path") > 0)
      IncludePaths = Variables["include-path"].as<std::vector<std::string>>();
    isl::ctx Context = allocateContextWithIncludePaths(IncludePaths);

    // run the cache model
    bullseyelib::run_model_new(Context,ConvertToProgramParameters(Variables));

  } catch (const boost::program_options::error &ex) {
    printf("-> exit(-1) option parsing error: %s\n", ex.what());
  }
  return 0;
}
