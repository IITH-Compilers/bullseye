// /******************************************************************************
//  *              libbullseye -  A library version of Bullseye *
//  ******************************************************************************/

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
#include "bullseyelib.h"

namespace bullseyelib {
bool check_path(std::string path) {
  std::ifstream f(path.c_str());
  return f.good();
}

std::map<int, std::string> compute_lines(std::string FileName,
                                         std::pair<long, long> ScopLoc) {
  std::map<int, std::string> Result;
  std::ifstream SourceFile;
  SourceFile.open(FileName);
  std::string Line;
  int LineNumber = 0;
  while (std::getline(SourceFile, Line)) {
    LineNumber++;
    if (SourceFile.tellg() > ScopLoc.first &&
        SourceFile.tellg() <= ScopLoc.second) {
      Result[LineNumber] = Line;
    }
  }
  SourceFile.close();
  return Result;
}

std::map<int, std::pair<long, long>> compute_offsets(std::string FileName) {
  std::map<int, std::pair<long, long>> Result;
  std::ifstream SourceFile;
  SourceFile.open(FileName);
  std::string Line;
  int LineNumber = 0;
  long Start = SourceFile.tellg();
  long End = 0;
  while (std::getline(SourceFile, Line)) {
    End = SourceFile.tellg();
    Result[++LineNumber] = std::make_pair(Start, End);
    Start = End;
  }
  SourceFile.close();
  return Result;
}

void print_scop(std::map<int, std::string> &Lines, int Start, int Stop) {
  // compute number of necessary digits
  int Width = std::to_string(Lines.end()->first).length();
  for (int i = Start; i < Stop; ++i) {
    std::cout << std::setw(Width) << std::right << std::to_string(i);
    std::cout << " " << Lines[i] << std::endl;
  }
}

std::vector<NamedLong>
ParametersFromDefineParametersIfPossible(ProgramParameters PP) {
  std::vector<NamedLong> Parameters;
  if (PP.DefineParameters.size() != 0) {
    printf("-> parsing the parameters...\n");
    std::vector<std::string> ParamStrings = PP.DefineParameters;
    for (auto ParamString : ParamStrings) {
      std::string Name;
      long Value;
      try {
        std::string Delimiter = "=";
        auto Split = ParamString.find(Delimiter);
        if (Split == std::string::npos) {
          throw std::runtime_error("did not find delimiter");
        }
        std::string Name = ParamString.substr(0, ParamString.find(Delimiter));
        long Value =
            std::stol(ParamString.substr(ParamString.find(Delimiter) + 1));
        printf("   - %s = %ld\n", Name.c_str(), Value);
        Parameters.push_back(std::make_pair(Name, Value));
      } catch (std::exception) {
        printf("-> exit(-1) failed to parse %s\n", ParamString.c_str());
        exit(-1);
      }
    }
    printf("-> done\n");
  }
  return Parameters;
}

void printCacheMissResults(ProgramParameters PP, CacheMissResults CMR) {
  // open the input file and seek the start of the scop
  std::map<int, std::string> Lines =
      compute_lines(PP.InputFile, CMR.Model.getScopLoc());
  std::map<int, std::pair<long, long>> Offsets = compute_offsets(PP.InputFile);
  long Position = Lines.begin()->first;
  std::string LineStart;
  std::string SingleLine;
  std::string DoubleLine;
  LineStart.resize(std::to_string(Lines.rbegin()->first).length() + 1, ' ');
  SingleLine.resize(80 - LineStart.length(), '-');
  DoubleLine.resize(80, '=');
  // print the access infos sorted by position
  size_t RefWidth = 16;
  std::map<long, std::vector<access_info>> AccessInfosByLn;
  std::map<std::string, access_info> AccessInfoByName;
  for (auto AccessInfos : CMR.Model.getAccessInfos()) {
    if (AccessInfos.second.empty())
      continue;
    AccessInfosByLn[AccessInfos.second[0].Line] = AccessInfos.second;
    for (auto AccessInfo : AccessInfos.second) {
      AccessInfoByName[AccessInfo.Name] = AccessInfo;
      RefWidth = std::max(RefWidth, AccessInfo.Name.length() + 1);
    }
  }
  // print the cache info access by access
  std::cout << DoubleLine << std::endl;
  std::cout << "                  relative number of cache misses (statement) "
            << std::endl;
  std::cout << DoubleLine << std::endl;
  for (auto AccessInfos : AccessInfosByLn) {
    // determine the last line of multiline
    int Next = AccessInfos.first;
    while (Offsets[Next].second < AccessInfos.second[0].Stop) {
      Next++;
    }
    // print the sources
    print_scop(Lines, Position, Next + 1);
    Position = Next + 1;
    // print header
    std::cout << LineStart << SingleLine << std::endl;
    std::cout << std::setw(RefWidth) << std::right << "ref";
    std::cout << "  ";
    std::cout << std::setw(6) << std::left << "type";
    std::cout << std::setw(10) << std::left << "comp[%]";
    for (int i = 1; i <= CMR.MachineModel.CacheSizes.size(); ++i) {
      std::string Capacity = "L" + std::to_string(i) + "[%]";
      std::cout << std::setw(10) << std::left << Capacity;
    }
    std::cout << std::setw(10) << std::left << "tot[%]";
    std::cout << std::setw(10) << std::left << "reuse[ln]";
    std::cout << std::endl;
    // print the accesses
    for (auto AccessInfo : AccessInfos.second) {
      // find the actual cache miss info
      auto Iter = std::find_if(
          CMR.CacheMisses.begin(), CMR.CacheMisses.end(),
          [&](NamedMisses Misses) { return Misses.first == AccessInfo.Name; });
      assert(Iter != CMR.CacheMisses.end());
      auto Compulsory = Iter->second.CompulsoryMisses;
      auto Capacity = Iter->second.CapacityMisses;
      auto Total = Iter->second.Total;
      // print the access info
      std::string Name = AccessInfo.Access;
      if (Name.length() > RefWidth)
        Name = Name.substr(0, RefWidth);
      std::cout << std::setw(RefWidth) << std::right << Name;
      std::cout << "  ";
      std::cout << std::setw(6) << std::left
                << (AccessInfo.ReadOrWrite == Read ? "rd" : "wr");
      std::cout << std::setw(10) << std::left << std::setprecision(5)
                << std::fixed
                << 100.0 * (double)Compulsory / (double)CMR.TotalAccesses;
      for (int i = 0; i < CMR.MachineModel.CacheSizes.size(); ++i) {
        std::cout << std::setw(10) << std::left << std::setprecision(5)
                  << std::fixed
                  << 100.0 * (double)Capacity[i] / (double)CMR.TotalAccesses;
      }
      std::cout << std::setw(10) << std::left << std::setprecision(5)
                << std::fixed
                << 100.0 * (double)Total / (double)CMR.TotalAccesses;
      // compute the reuse line numbers
      auto Conflicts = CMR.Model.getConflicts()[AccessInfo.Name];
      // compute the reuse line numbers
      std::vector<int> ReuseLines;
      for (auto Conflict : Conflicts) {
        ReuseLines.push_back(AccessInfoByName[Conflict].Line);
      }
      // sort the line numbers and remove duplicates
      std::sort(ReuseLines.begin(), ReuseLines.end());
      auto Last = std::unique(ReuseLines.begin(), ReuseLines.end());
      for (auto Iter = ReuseLines.begin(); Iter != Last;) {
        std::cout << *Iter;
        if (++Iter != Last)
          std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << LineStart << SingleLine << std::endl;
  }
  print_scop(Lines, Position, Lines.rbegin()->first + 1);
  // print the scop info
  std::cout << DoubleLine << std::endl;
  std::cout << "                     absolute number of cache misses (SCOP)"
            << std::endl;
  std::cout << DoubleLine << std::endl;
  std::cout.imbue(std::locale(""));
  std::cout << std::setw(16) << std::left << "compulsory:";
  std::cout << std::setw(20) << std::right << CMR.TotalCompulsory << std::endl;
  for (int i = 1; i <= CMR.MachineModel.CacheSizes.size(); ++i) {
    std::string Capacity = "capacity (L" + std::to_string(i) + "):";
    std::cout << std::setw(16) << std::left << Capacity;
    std::cout << std::setw(20) << std::right << CMR.TotalCapacity[i - 1]
              << std::endl;
  }
  std::cout << std::setw(16) << std::left << "total:";
  std::cout << std::setw(20) << std::right << CMR.TotalAccesses << std::endl;
  std::cout << std::setw(16) << std::left << "flops:";
  std::cout << std::setw(20) << std::right << CMR.TotalFlopCount_ << std::endl;
  std::cout << DoubleLine << std::endl;
}

CacheMissResults run_model_new(isl::ctx Context, ProgramParameters PP) {
  // allocate the machine model with default values
  machine_model MachineModel = {PP.LineSize, PP.CacheSizes};
  model_options ModelOptions = {PP.ComputeBounds};
  printf("-> setting up cache levels\n");
  std::sort(MachineModel.CacheSizes.begin(), MachineModel.CacheSizes.end());
  for (auto CacheSize : MachineModel.CacheSizes) {
    if (CacheSize % 1024 == 0) {
      printf("   - %ldkB with %ldB cache lines\n", CacheSize / 1024,
             MachineModel.CacheLineSize);
    } else {
      printf("   - %ldB with %ldB cache lines\n", CacheSize,
             MachineModel.CacheLineSize);
    }
  }
  printf("-> done\n");
  // compute the total time
  auto StartExecution = std::chrono::high_resolution_clock::now();
  // allocate the cache model and compile the program
  HayStack Model(Context, MachineModel, ModelOptions);
  if (PP.ScopFunction.empty()) {
    Model.compileProgram(PP.InputFile);
  } else {
    Model.compileProgram(PP.InputFile, PP.ScopFunction);
  }
  // parsing the parameters
  std::vector<NamedLong> Parameters =
      ParametersFromDefineParametersIfPossible(PP);
  // run the preprocessing
  printf("-> start processing...\n");
  auto Start = std::chrono::high_resolution_clock::now();
  auto Start_Model = std::chrono::high_resolution_clock::now();
  Model.initModel(Parameters);
  auto Stop_Model = std::chrono::high_resolution_clock::now();
  double TotalEvaluation_init =
      std::chrono::duration<double, std::milli>(Stop_Model - Start_Model)
          .count();
  printf("-> Intialization done after (%.2fms)\n", TotalEvaluation_init);
  // execute the cache model
  auto CacheMisses = Model.countCacheMisses();
  auto Stop = std::chrono::high_resolution_clock::now();
  double TotalEvaluation =
      std::chrono::duration<double, std::milli>(Stop - Start).count();
  printf("-> done after (%.2fms)\n", TotalEvaluation);
  // collect and print result
  long TotalAccesses = 0;
  long TotalCompulsory = 0;
  std::vector<long> TotalCapacity(MachineModel.CacheSizes.size(), 0);
  // sum the cache misses for all accesses
  for (auto &CacheMiss : CacheMisses) {
    TotalAccesses += CacheMiss.second.Total;
    TotalCompulsory += CacheMiss.second.CompulsoryMisses;
    std::transform(TotalCapacity.begin(), TotalCapacity.end(),
                   CacheMiss.second.CapacityMisses.begin(),
                   TotalCapacity.begin(), std::plus<long>());
  };
  CacheMissResults CMR(TotalCompulsory, TotalAccesses, TotalCapacity,
                       MachineModel, ModelOptions, Model, Model.countFlops(), CacheMisses);
  printCacheMissResults(PP, CMR);
  return CMR;
}

} // namespace bullseyelib
