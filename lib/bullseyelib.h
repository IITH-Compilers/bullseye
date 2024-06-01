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
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

#include <isl/options.h>

#include "Definitions.h"
#include "HayStack.h"
#include "Timer.h"

namespace bullseyelib {
const int CACHE_LINE_SIZE = 64;

// default
const int CACHE_SIZE1 = 32 * 1024;
const int CACHE_SIZE2 = 512 * 1024;

struct ProgramParameters {
  std::vector<long> CacheSizes;
  long LineSize;
  std::string InputFile;
  std::vector<std::string> IncludePath;
  std::vector<std::string> DefineParameters;
  std::string ScopFunction;
  bool ComputeBounds; 

  ProgramParameters() = delete;

  ProgramParameters(std::string inputfile,
                    std::vector<long> cachesizes = {CACHE_SIZE1, CACHE_SIZE2},
                    long linesize = CACHE_LINE_SIZE,
                    std::vector<std::string> includepath = {},
                    std::vector<std::string> defineparameters = {},
                    std::string scopfunction = "", bool computebounds = false)
      : CacheSizes(cachesizes), LineSize(linesize), InputFile(inputfile),
        IncludePath(includepath), DefineParameters(defineparameters),
        ComputeBounds(computebounds) {}
};

struct CacheMissResults {
  long TotalCompulsory;
  long TotalAccesses;
  std::vector<long> TotalCapacity;
  machine_model MachineModel;
  model_options ModelOptions;
  HayStack Model;
  long long TotalFlopCount_;
  std::vector<NamedMisses> CacheMisses;

  CacheMissResults(long totalcompulsory, long totalaccesses,
                   std::vector<long> totalcapacity, machine_model mm,
                   model_options mo, HayStack model, long long totalflopcount,
                   std::vector<NamedMisses> cm)
      : TotalCompulsory(totalcompulsory), TotalAccesses(totalaccesses),
        TotalCapacity(totalcapacity), MachineModel(mm), ModelOptions(mo),
        Model(model), TotalFlopCount_(totalflopcount) ,CacheMisses(cm) {}
};

std::vector<NamedLong>
ParametersFromDefineParametersIfPossible(ProgramParameters PP);

bool check_path(std::string path);

std::map<int, std::string> compute_lines(std::string FileName,
                                         std::pair<long, long> ScopLoc);

std::map<int, std::pair<long, long>> compute_offsets(std::string FileName);

void print_scop(std::map<int, std::string> &Lines, int Start, int Stop);

CacheMissResults run_model_new(isl::ctx Context, ProgramParameters PP);

std::vector<NamedLong>
ParametersFromDefineParametersIfPossible(ProgramParameters PP);

void printCacheMissResults(ProgramParameters PP, CacheMissResults CMR);

} // namespace bullseyelib

#endif
