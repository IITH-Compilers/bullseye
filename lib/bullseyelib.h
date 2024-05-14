/******************************************************************************
 *              libbullseye -  A library version of Bullseye                  *
 ******************************************************************************/
/*
 * Copyright (c) 2022, BullsEye
   Adopted from HayStack (PLDI 2019)
*/

/*
 * Copyright (c) 2022, BullsEye
*/

#ifndef _BULLSEYELIB_H_
#define _BULLSEYELIB_H_

// #include<iostream>

// void bullseye();


#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

#include <isl/options.h>

#include "HayStack.h"
#include "Timer.h"

namespace po = boost::program_options;

namespace bullseyelib {
  bool check_path(std::string path);

  std::map<int, std::string> compute_lines(std::string FileName, std::pair<long, long> ScopLoc);

  std::map<int, std::pair<long, long>> compute_offsets(std::string FileName);

  void print_scop(std::map<int, std::string> &Lines, int Start, int Stop);

  void run_model(isl::ctx Context, po::variables_map Variables);
}

#endif