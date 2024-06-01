/*
 * Copyright (c) 2022, BullsEye
 */

#ifndef _BULLSEYE_H_
#define _BULLSEYE_H_

#include <boost/math/common_factor.hpp>
#include <cmath>
#include <ilcplex/ilocplex.h>
#include <isl/isl-noexceptions.h>
#include <isl/vertices.h>
#include <map>
#include <string>
#include <vector>

#include "Definitions.h"

// Global variables
struct oct {
  std::map<long, long> coeff;
  long constant;
};
struct data_poly {
  int index1;
  int index2;
  int dims;
  size_t psize;
  std::vector<std::vector<long>> transpose;
  int limit;
  std::vector<std::vector<long>> hmatrix;
  isl::local_space ls;
};

/*
    Implementation of BullsEye
*/
class BullsEye {
public:
  BullsEye() = delete;
  BullsEye(const BullsEye &other) = default;

  // additonal functions for preprocessing and computations
  static std::vector<isl::qpolynomial>
  _combinations_(isl_pw_qpolynomial *poly, int N, int K,
                 std::vector<isl::qpolynomial> mon);
  static long _factorial_(long number);
  static long _combineNCR_(long n, long r);
  static isl::pw_qpolynomial clean_poly_lcm(isl::qpolynomial qcopy);

  // get the LP bounds for the given indices and boundtype
  static struct oct
  findLPIntervalBounds(int index1, int index2, int boundtype1, int boundtype2,
                       int varsize, int psize,
                       std::vector<std::vector<long>> &transpose, long limit,
                       std::vector<std::vector<long>> &hmatrix);
  static struct oct findLPBounds(int index1, int index2, int boundtype1,
                                 int boundtype2, int varsize, int psize,
                                 std::vector<std::vector<long>> &transpose,
                                 long limit,
                                 std::vector<std::vector<long>> &hmatrix);
  static struct oct findLPBoundsNeg(int index1, int index2, int boundtype1,
                                    int boundtype2, int varsize, int psize,
                                    std::vector<std::vector<long>> &transpose,
                                    long limit,
                                    std::vector<std::vector<long>> &hmatrix);

  // Methods for Vertex based Counting
  static isl_stat vertex_count(__isl_take isl_vertex *vertex, void *user);
  static isl_stat vertex_count_interval(__isl_take isl_vertex *vertex,
                                        void *user);
  static isl_stat vertex_test(__isl_take isl_vertex *vertex, void *user);

  // Extract Affine expression
  static isl::pw_aff extractAffineExpression(piece Piece);
  static std::vector<long> countAffineDimensions(piece Piece,
                                                 std::vector<long> Limits);
  static isl::aff extractAffineExpression(isl::qpolynomial Polynomial,
                                          isl::set Domain,
                                          std::map<int, long> Values);

  // Method to set the HMatrix values
  static void setcolumn(std::string bmask, int column, long long coeff,
                        std::vector<std::vector<long>> &hm);

  // compute approximate capacity misses
  static std::vector<long>
  calculateApproximateCapacityMisses(piece Piece, std::vector<int> NonAffine,
                                     std::vector<int> Affine,
                                     std::vector<long> cache_size);
};

#endif