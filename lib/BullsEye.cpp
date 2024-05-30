/*
 * Copyright (c) 2022, BullsEye
 */

#include <algorithm>
#include <cassert>
#include <limits>

#include "Access.h"
#include "BullsEye.h"
#include "Timer.h"
#include "barvinok/isl.h"
#include "isl-helpers.h"

struct check_poly {
  isl::pw_qpolynomial _polynomial;
  isl::basic_set input_domain;
  long cache_size;
};

std::vector<isl::qpolynomial>
BullsEye::_combinations_(isl_pw_qpolynomial *poly, int N, int K,
                         std::vector<isl::qpolynomial> mon) {

  std::vector<isl::qpolynomial> temp(mon);
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);
  isl::space sp =
      isl::manage(isl_space_copy(isl_pw_qpolynomial_get_domain_space(poly)));
  isl::local_space ls = isl::manage(isl_local_space_from_space(sp.release()));
  isl::aff a1 = isl::manage(isl_aff_zero_on_domain(ls.release()));
  a1 = a1.set_constant_si(1);
  isl::qpolynomial qtemp = isl::manage(isl_qpolynomial_from_aff(a1.release()));

  temp.push_back(qtemp);
  do {
    isl::qpolynomial qtemp2(qtemp);
    for (int i = 0; i < N; ++i) {
      if (bitmask[i]) {
        isl::qpolynomial qtemp3(mon[i]);
        qtemp2 = isl::manage(
            isl_qpolynomial_mul(qtemp2.release(), qtemp3.release()));
      }
    }
    temp.push_back(qtemp2);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  for (auto i : mon) {
    isl::qpolynomial qtemp3(i);
    isl::qpolynomial qtemp4(i);
    qtemp3 =
        isl::manage(isl_qpolynomial_mul(qtemp4.release(), qtemp3.release()));
    temp.push_back(qtemp3);
  }
  return temp;
}

// Return set of Interval constraints
struct oct BullsEye::findLPIntervalBounds(
    int index1, int index2, int boundtype1, int boundtype2, int varsize,
    int psize, std::vector<std::vector<long>> &transpose, long limit,
    std::vector<std::vector<long>> &hmatrix) {

  struct oct constraint;
  constraint.constant = 0;
  IloEnv env;
  IloModel model(env);
  IloFloatVarArray x(env, psize - 1, 0, IloInfinity, ILOFLOAT);

  // Objective function
  IloExpr exp0(env);
  int sizetrans = transpose.size();
  for (int s = 0; s < sizetrans - 1; s++) {
    int temp = 0;
    temp += transpose[s][0];
    temp += boundtype1 * transpose[s][index1 + 1];
    temp += boundtype2 * transpose[s][index2 + 1];
    exp0 += x[s] * temp;
  }
  exp0 += transpose[sizetrans - 1][0];
  exp0 += boundtype1 * transpose[sizetrans - 1][index1 + 1];
  exp0 += boundtype2 * transpose[sizetrans - 1][index2 + 1];
  model.add(IloMinimize(env, exp0));
  IloExpr cons0(env, 0);
  IloIntExpr cons1(env, 1);
  IloIntExpr consNeg(env, -1);

  // Add Constraints
  // First constraints to cancel the non-affine terms
  for (int s = varsize + 1; s < hmatrix.size(); s++) {
    IloExpr exp1(env);
    for (int d = 0; d < hmatrix[s].size() - 1; d++) {
      exp1 += x[d] * hmatrix[s][d];
    }
    exp1 += hmatrix[s][hmatrix[s].size() - 1];
    model.add(exp1 == cons0);
  }

  // Add constraints for Octagons
  IloExpr exp2(env);
  for (int d = 0; d < hmatrix[index1 + 1].size() - 1; d++) {
    exp2 += x[d] * hmatrix[index1 + 1][d];
  }
  exp2 += hmatrix[index1 + 1][hmatrix[index1 + 1].size() - 1];

  IloExpr exp3(env);
  for (int d = 0; d < hmatrix[index2 + 1].size() - 1; d++) {
    exp3 += x[d] * hmatrix[index2 + 1][d];
  }
  exp3 += hmatrix[index2 + 1][hmatrix[index2 + 1].size() - 1];
  model.add(exp3 == cons0);

  for (int i = 0; i < psize - 1; i++) {
    IloExpr exp2(env);
    exp2 += x[i];
    model.add(exp2 >= cons0);
  }

  IloCplex cplex(model);
  cplex.setOut(env.getNullStream());
  cplex.extract(model);

  cplex.solve();
  exp0.end();
  exp2.end();
  exp3.end();

  if (cplex.getStatus() == IloAlgorithm::Infeasible) {
    return constraint;
  }
  double bounds = 0;

  if (cplex.getStatus() == IloAlgorithm::Optimal) {
    IloFloatArray vals(env);
    cplex.getValues(vals, x);
    for (int i = 1; i <= varsize; i++) {
      double const1 = 0;
      for (int d = 0; d < hmatrix[i].size() - 1; d++) {
        const1 += vals[d] * hmatrix[i][d];
      }
      const1 += hmatrix[i][hmatrix[i].size() - 1];
      constraint.coeff[i] = const1;
    }
    for (int d = 0; d < hmatrix[0].size() - 1; d++) {
      bounds += vals[d] * hmatrix[0][d];
    }
    bounds += hmatrix[0][hmatrix[0].size() - 1];
    constraint.constant = ceil(bounds);
  }
  env.end();
  return constraint;
}

// Return set of octagonal constraints for Positive bounds
struct oct BullsEye::findLPBounds(int index1, int index2, int boundtype1,
                                  int boundtype2, int varsize, int psize,
                                  std::vector<std::vector<long>> &transpose,
                                  long limit,
                                  std::vector<std::vector<long>> &hmatrix) {
  struct oct constraint;
  constraint.constant = 0;
  IloEnv env;
  IloModel model(env);
  IloFloatVarArray x(env, psize - 1, 0, IloInfinity, ILOFLOAT);

  // Objective function
  IloExpr exp0(env);
  int sizetrans = transpose.size();
  for (int s = 0; s < sizetrans - 1; s++) {
    int temp = 0;
    temp += transpose[s][0];
    temp += boundtype1 * transpose[s][index1 + 1];
    temp += boundtype2 * transpose[s][index2 + 1];
    exp0 += x[s] * temp;
  }
  exp0 += transpose[sizetrans - 1][0];
  exp0 += boundtype1 * transpose[sizetrans - 1][index1 + 1];
  exp0 += boundtype2 * transpose[sizetrans - 1][index2 + 1];
  model.add(IloMinimize(env, exp0));
  IloExpr cons0(env, 0);
  IloIntExpr cons1(env, 1);
  IloIntExpr consNeg(env, -1);

  // Add Constraints
  // First constraints to cancel the non-affine terms
  for (int s = varsize + 1; s < hmatrix.size(); s++) {
    IloExpr exp1(env);
    for (int d = 0; d < hmatrix[s].size() - 1; d++) {
      exp1 += x[d] * hmatrix[s][d];
    }
    exp1 += hmatrix[s][hmatrix[s].size() - 1];
    model.add(exp1 == cons0);
  }

  // Add constraints for Octagons
  IloExpr exp2(env);
  for (int d = 0; d < hmatrix[index1 + 1].size() - 1; d++) {
    exp2 += x[d] * hmatrix[index1 + 1][d];
  }
  exp2 += hmatrix[index1 + 1][hmatrix[index1 + 1].size() - 1];

  IloExpr exp3(env);
  for (int d = 0; d < hmatrix[index2 + 1].size() - 1; d++) {
    exp3 += x[d] * hmatrix[index2 + 1][d];
  }
  exp3 += hmatrix[index2 + 1][hmatrix[index2 + 1].size() - 1];
  model.add(exp3 == exp2);

  for (int i = 0; i < psize - 1; i++) {
    IloExpr exp2(env);
    exp2 += x[i];
    model.add(exp2 >= cons0);
  }

  IloCplex cplex(model);
  cplex.setOut(env.getNullStream());
  cplex.extract(model);

  cplex.solve();
  exp0.end();
  exp2.end();
  exp3.end();

  if (cplex.getStatus() == IloAlgorithm::Infeasible) {
    return constraint;
  }
  double bounds = 0;

  if (cplex.getStatus() == IloAlgorithm::Optimal) {
    IloFloatArray vals(env);
    cplex.getValues(vals, x);
    for (int i = 1; i <= varsize; i++) {
      double const1 = 0;
      for (int d = 0; d < hmatrix[i].size() - 1; d++) {
        const1 += vals[d] * hmatrix[i][d];
      }
      const1 += hmatrix[i][hmatrix[i].size() - 1];
      constraint.coeff[i] = const1;
    }
    for (int d = 0; d < hmatrix[0].size() - 1; d++) {
      bounds += vals[d] * hmatrix[0][d];
    }
    bounds += hmatrix[0][hmatrix[0].size() - 1];
    constraint.constant = ceil(bounds);
  }
  env.end();
  return constraint;
}

// Return set of octagonal constraints for negative bounds
struct oct BullsEye::findLPBoundsNeg(int index1, int index2, int boundtype1,
                                     int boundtype2, int varsize, int psize,
                                     std::vector<std::vector<long>> &transpose,
                                     long limit,
                                     std::vector<std::vector<long>> &hmatrix) {
  struct oct constraint;
  constraint.constant = 0;
  IloEnv env;
  IloModel model(env);
  IloFloatVarArray x(env, psize - 1, 0, IloInfinity, ILOFLOAT);

  // Objective function
  IloExpr exp0(env);
  int sizetrans = transpose.size();
  for (int s = 0; s < sizetrans - 1; s++) {
    int temp = 0;
    temp += transpose[s][0];
    temp += boundtype1 * transpose[s][index1 + 1];
    temp += boundtype2 * transpose[s][index2 + 1];
    exp0 += x[s] * temp;
  }
  exp0 += transpose[sizetrans - 1][0];
  exp0 += boundtype1 * transpose[sizetrans - 1][index1 + 1];
  exp0 += boundtype2 * transpose[sizetrans - 1][index2 + 1];
  model.add(IloMinimize(env, exp0));
  IloExpr cons0(env, 0);
  IloIntExpr cons1(env, 1);
  IloIntExpr consNeg(env, -1);

  // Add Constraints
  // First constraints to cancel the non-affine terms
  for (int s = varsize + 1; s < hmatrix.size(); s++) {
    IloExpr exp1(env);
    for (int d = 0; d < hmatrix[s].size() - 1; d++) {
      exp1 += x[d] * hmatrix[s][d];
    }
    exp1 += hmatrix[s][hmatrix[s].size() - 1];
    model.add(exp1 == cons0);
  }

  // Add constraints for Octagons
  IloExpr exp2(env);
  for (int d = 0; d < hmatrix[index1 + 1].size() - 1; d++) {
    exp2 += x[d] * hmatrix[index1 + 1][d];
  }
  exp2 += hmatrix[index1 + 1][hmatrix[index1 + 1].size() - 1];

  IloExpr exp3(env);
  for (int d = 0; d < hmatrix[index2 + 1].size() - 1; d++) {
    exp3 += x[d] * hmatrix[index2 + 1][d];
  }
  exp3 += hmatrix[index2 + 1][hmatrix[index2 + 1].size() - 1];
  model.add(exp2 == -exp3);

  for (int i = 0; i < psize - 1; i++) {
    IloExpr exp2(env);
    exp2 += x[i];
    model.add(exp2 >= cons0);
  }

  IloCplex cplex(model);
  cplex.setOut(env.getNullStream());
  cplex.extract(model);

  cplex.solve();

  if (cplex.getStatus() == IloAlgorithm::Infeasible) {
    return constraint;
  }
  double bounds = 0;

  if (cplex.getStatus() == IloAlgorithm::Optimal) {
    IloFloatArray vals(env);
    cplex.getValues(vals, x);
    for (int i = 1; i <= varsize; i++) {
      double const1 = 0;
      for (int d = 0; d < hmatrix[i].size() - 1; d++) {
        const1 += vals[d] * hmatrix[i][d];
      }
      const1 += hmatrix[i][hmatrix[i].size() - 1];
      constraint.coeff[i] = const1;
    }
    for (int d = 0; d < hmatrix[0].size() - 1; d++) {
      bounds += vals[d] * hmatrix[0][d];
    }
    bounds += hmatrix[0][hmatrix[0].size() - 1];
    constraint.constant = ceil(bounds);
  }
  env.end();
  return constraint;
}

std::vector<isl::constraint> constraint_list;

// Return factorial of a integer
long BullsEye::_factorial_(long n) {
  int res = 1;
  for (int i = 2; i <= n; i++)
    res = res * i;
  return res;
}

// Return nCr value for input n and r
long BullsEye::_combineNCR_(long n, long r) {
  return _factorial_(n) / (_factorial_(r) * _factorial_(n - r));
}

isl::pw_qpolynomial BullsEye::clean_poly_lcm(isl::qpolynomial qcopy) {

  isl::space sp = isl::manage(
      isl_space_copy(isl_qpolynomial_get_domain_space(qcopy.get())));
  isl::pw_qpolynomial pwq1 = isl::manage(isl_pw_qpolynomial_from_qpolynomial(
      isl_qpolynomial_zero_on_domain(sp.get())));
  isl::local_space ls = isl::manage(isl_local_space_from_space(sp.release()));

  std::vector<long> Mulipliers;

  auto findMultipliers = [&](isl::term Term) {
    // get the divisor variable
    long multiplier = 1;
    long multip = 1;
    int max_exp = 1;
    for (int j = 0; j < Term.dim(isl::dim::set); j++) {
      if (max_exp < Term.get_exp(isl::dim::set, j)) {
        max_exp = Term.get_exp(isl::dim::set, j);
      }
    }
    isl::qpolynomial k_new = isl::qpolynomial::from_term(Term);
    isl::pw_qpolynomial k_pw =
        isl::manage(isl_pw_qpolynomial_from_qpolynomial(k_new.release()));
    std::string s_poly = k_pw.to_str();
    k_pw.release();

    size_t start_pos = s_poly.find('/');
    if (start_pos != std::string::npos) {
      std::string copy_str(s_poly);
      copy_str = copy_str.substr(start_pos + 1, copy_str.length() - 1);
      std::size_t found = copy_str.find_first_not_of("0123456789");
      char c = ' ';
      if (found != std::string::npos) {
        c = copy_str[found];
        copy_str = copy_str.substr(0, found);
        multiplier = std::stoi(copy_str);
        multip = multiplier;
        size_t pos = std::string::npos;
        if (!(c == ' ') && (pos = s_poly.find("^")) != std::string::npos) {
          max_exp = 2;
          multiplier = (int)(pow(multiplier, max_exp) + 0.5);
        }
      }

      Mulipliers.push_back(multiplier);

      if (c == ')') {
        std::string toErase = "floor";
        size_t pos = std::string::npos;

        while ((pos = s_poly.find(toErase)) != std::string::npos) {
          s_poly.erase(pos, toErase.length());
        }
        size_t start_pos2 = s_poly.find('/');
        std::string copy_str2(s_poly);

        copy_str2 = copy_str2.substr(start_pos2, copy_str2.length() - 1);

        std::string subc = copy_str2.substr(0, found + 1);
        pos = std::string::npos;
        int count_floors = 0;
        while ((pos = s_poly.find(subc)) != std::string::npos) {
          count_floors++;
          s_poly.erase(pos, subc.length());
        }
        std::string numerator = "1";
        if (max_exp > 1) {
          subc = "/" + std::to_string(multiplier);
        } else if (count_floors > 1) {
          multip = (int)(pow(multip, count_floors) + 0.5);
          subc = "/" + std::to_string(multip);
        }
        subc = numerator.append(subc);
        size_t start_pos3 = s_poly.find('>');
        if (start_pos3 != std::string::npos) {
          std::string mult = subc + "*";
          s_poly.insert(start_pos3 + 1, std::string(mult));
        }
      }
    }
    isl::pw_qpolynomial pwq2 = isl::manage(isl_pw_qpolynomial_read_from_str(
        qcopy.get_ctx().get(), s_poly.c_str()));

    pwq1 = isl::manage(isl_pw_qpolynomial_add(pwq1.release(), pwq2.release()));

    return isl::stat::ok();
  };
  qcopy.foreach_term(findMultipliers);

  long lcm = 1;
  for (auto mul : Mulipliers) {
    lcm = boost::math::lcm(mul, lcm);
  }

  isl::aff a1 = isl::manage(isl_aff_zero_on_domain(ls.release()));
  a1 = a1.set_constant_si(lcm);
  isl::qpolynomial qtemp = isl::manage(isl_qpolynomial_from_aff(a1.release()));
  pwq1 = pwq1.mul(
      isl::manage(isl_pw_qpolynomial_from_qpolynomial(qtemp.release())));
  qcopy.release();
  return pwq1;
}

isl_stat BullsEye::vertex_count(__isl_take isl_vertex *vertex, void *user) {

  constraint_list.clear();
  struct data_poly *data_p = (struct data_poly *)user;
  isl_multi_val *mval =
      isl_multi_aff_get_constant_multi_val(isl_vertex_get_expr(vertex));

  isl::multi_val mv(isl::manage(mval));
  std::map<long, long long> Values;
  auto index1 = data_p->index1;
  auto index2 = data_p->index2;
  auto dims = data_p->dims;
  auto psize = data_p->psize;
  auto transpose = data_p->transpose;
  auto hmatrix = data_p->hmatrix;
  auto lspace = data_p->ls;

  for (int i = 0; i < dims; ++i) {
    auto Value = mv.get_val(i);
    Values[i] = Value.get_num_si();
  }
  struct oct b1 = findLPBounds(index1, index2, Values[index1], Values[index2],
                               dims, psize, transpose, 512, hmatrix);
  struct oct b2 =
      findLPBoundsNeg(index1, index2, Values[index1], Values[index2], dims,
                      psize, transpose, 512, hmatrix);

  if (b1.coeff[index1 + 1] != 0 && b1.coeff[index2 + 1] != 0) {

    isl::constraint lowerConsx = isl::constraint::alloc_inequality(lspace);
    lowerConsx = lowerConsx.set_coefficient_si(isl::dim::set, index1,
                                               b1.coeff[index1 + 1]);
    lowerConsx = lowerConsx.set_coefficient_si(isl::dim::set, index2,
                                               b1.coeff[index2 + 1]);
    isl::val v = isl::manage(
        isl_val_int_from_si(lowerConsx.get_ctx().get(), b1.constant));
    lowerConsx = lowerConsx.set_constant_val(v);
    constraint_list.push_back(lowerConsx);
  }

  if (b2.coeff[index1 + 1] != 0 && b2.coeff[index2 + 1] != 0) {

    isl::constraint sec = isl::constraint::alloc_inequality(lspace);
    sec = sec.set_coefficient_si(isl::dim::set, index1, b2.coeff[index1 + 1]);
    sec = sec.set_coefficient_si(isl::dim::set, index2, b2.coeff[index2 + 1]);
    isl::val v =
        isl::manage(isl_val_int_from_si(sec.get_ctx().get(), b2.constant));
    sec = sec.set_constant_val(v);
    constraint_list.push_back(sec);
  }

  isl_vertex_free(vertex);
  return isl_stat_ok;
}

isl_stat BullsEye::vertex_count_interval(__isl_take isl_vertex *vertex,
                                         void *user) {

  constraint_list.clear();
  struct data_poly *data_p = (struct data_poly *)user;
  isl_multi_val *mval =
      isl_multi_aff_get_constant_multi_val(isl_vertex_get_expr(vertex));

  isl::multi_val mv(isl::manage(mval));
  std::map<long, long long> Values;
  auto index1 = data_p->index1;
  auto index2 = data_p->index2;
  auto dims = data_p->dims;
  auto psize = data_p->psize;
  auto transpose = data_p->transpose;
  auto hmatrix = data_p->hmatrix;
  auto lspace = data_p->ls;

  for (int i = 0; i < dims; ++i) {
    auto Value = mv.get_val(i);
    Values[i] = Value.get_num_si();
  }

  struct oct b1 =
      findLPIntervalBounds(index1, index2, Values[index1], Values[index2], dims,
                           psize, transpose, 512, hmatrix);

  if (b1.coeff[index1 + 1] != 0) {

    isl::constraint lowerConsx = isl::constraint::alloc_inequality(lspace);
    lowerConsx = lowerConsx.set_coefficient_si(isl::dim::set, index1,
                                               b1.coeff[index1 + 1]);
    lowerConsx = lowerConsx.set_coefficient_si(isl::dim::set, index2,
                                               b1.coeff[index2 + 1]);
    isl::val v = isl::manage(
        isl_val_int_from_si(lowerConsx.get_ctx().get(), b1.constant));
    lowerConsx = lowerConsx.set_constant_val(v);

    constraint_list.push_back(lowerConsx);
  }

  isl_vertex_free(vertex);
  return isl_stat_ok;
}

int ver_num = 0;
std::vector<bool> count_poly(10, 0);

// Test a input vertex for a cache miss
isl_stat BullsEye::vertex_test(__isl_take isl_vertex *vertex, void *user) {

  struct check_poly *data_p = (struct check_poly *)user;
  isl_multi_val *mval =
      isl_multi_aff_get_constant_multi_val(isl_vertex_get_expr(vertex));

  auto poly = data_p->_polynomial;
  auto domain = data_p->input_domain;
  auto cache_size = data_p->cache_size;
  isl::multi_val mv(isl::manage(mval));
  isl::point Point = isl::manage(isl_point_zero(isl_multi_val_get_space(mval)));

  for (int i = 0; i < Point.get_space().dim(isl::dim::set); ++i) {
    auto Value = mv.get_val(i);
    if (!Value.is_int())
      Value = Value.floor();
    Point = Point.set_coordinate_val(isl::dim::set, i, Value);
    Value.release();
  }

  isl::val StackDistance =
      isl::manage(isl_pw_qpolynomial_eval(poly.release(), Point.release()));

  if (isl::get_value(StackDistance) > cache_size)
    count_poly[ver_num] = true;
  else
    count_poly[ver_num] = false;
  ver_num++;

  mv.release();
  StackDistance.release();
  domain.release();
  isl_vertex_free(vertex);

  return isl_stat_ok;
}

// Count Approximate cache misses using Handelman-Octagon counting of stack
// distance polynomials
std::vector<long> BullsEye::calculateApproximateCapacityMisses(
    piece Piece, std::vector<int> NonAffine, std::vector<int> Affine,
    std::vector<long> cache_size) {

  isl::pw_qpolynomial _polynomial =
      isl::manage(isl_pw_qpolynomial_from_qpolynomial(
          isl_qpolynomial_copy(Piece.Polynomial.get())));
  ;
  isl::set input_domain = isl::manage(isl_set_copy(Piece.Domain.get()));
  isl::pw_aff expression;
  std::vector<long> Misses;

  bool interval_bounds = false;
  isl_set *trial_set_new = input_domain.release();
  isl_ctx *ctx = isl_set_get_ctx(trial_set_new);

  trial_set_new = isl_set_project_out(trial_set_new, isl_dim_set,
                                      isl_set_n_dim(trial_set_new) - 1, 1);
  isl_set *trial_set = isl_set_remove_divs(isl_set_copy(trial_set_new));
  _polynomial =
      _polynomial.drop_dims(isl::dim::in, isl_set_n_dim(trial_set_new), 1);
  long long totalMissCount = 0;
  long long cache_misses = 0;
  int cache_number = 0;
  unsigned enum_length = 25;

  // For each cache level
  for (auto cache_limit : cache_size) {

    // Initialize cache misses per cache level
    cache_misses = 0;
    // For Domain of dimensions greater than three, apply Sparse enumeration
    if (isl_set_n_dim(trial_set) > 3) {
      auto countCacheMisses = [&](isl::basic_set All) {
        isl::set Enumeration = All;
        for (auto Iter = Affine.rbegin(); Iter != Affine.rend(); Iter++) {
          Enumeration = Enumeration.project_out(isl::dim::set, *Iter, 1);
        }
        int iterCount = 0;
        auto countSubdomain = [&](isl::point Point) {
          if (iterCount % enum_length == 0) {
            // compute the updated domain of the piece
            std::map<int, long> Values;
            auto Domain = All;
            for (int i = 0; i < Point.get_space().dim(isl::dim::set); ++i) {
              auto Value = Point.get_coordinate_val(isl::dim::set, i);
              Domain = Domain.fix_val(isl::dim::set, NonAffine[i], Value);
              Values[NonAffine[i]] = isl::get_value(Value);
            }
            auto Expression = BullsEye::extractAffineExpression(
                Piece.Polynomial, Piece.Domain, Values);
            if (Expression.is_cst()) {
              long Size = isl::cardinality(Domain);
              auto Constant = Expression.get_constant_val();
              auto Misses =
                  Constant.gt(isl::val(Constant.get_ctx(), cache_limit)) ? Size
                                                                         : 0;
              cache_misses += Misses * enum_length;
            } else {
              input_domain = Domain;
              expression = Expression;
              auto Misses = BullsEye::countAffineDimensions(Piece, cache_size);
              std::transform(
                  Misses.begin(), Misses.end(), Misses.begin(),
                  [&enum_length](auto &c) { return c * enum_length; });
              cache_misses += Misses[cache_number];
            }
          } else {
            iterCount++;
          }
          return isl::stat::ok();
        };

        Enumeration.foreach_point(countSubdomain);
        return isl::stat::ok();
      };
      input_domain = isl::manage(isl_set_make_disjoint(input_domain.release()));
      input_domain.foreach_basic_set(countCacheMisses);
      Misses.push_back(cache_misses);
      cache_number++;

    } else {

      totalMissCount = 0;
      isl::qpolynomial tqp, tempqp;
      isl::pw_qpolynomial qcopy = isl::manage(_polynomial.copy());
      isl::space sp = isl::manage(
          isl_space_copy(isl_pw_qpolynomial_get_domain_space(qcopy.copy())));
      isl::local_space ls =
          isl::manage(isl_local_space_from_space(sp.release()));
      isl::aff affine_limit = isl::manage(isl_aff_zero_on_domain(ls.release()));
      affine_limit = affine_limit.set_constant_si(cache_limit);
      isl::qpolynomial qtemp =
          isl::manage(isl_qpolynomial_from_aff(affine_limit.release()));
      isl::pw_qpolynomial pwq1 =
          isl::manage(isl_pw_qpolynomial_from_qpolynomial(qtemp.release()));
      qcopy = qcopy.sub(isl::manage(pwq1.release()));

      auto analyzePoly = [&](isl::set s, isl::qpolynomial qp) {
        tqp = qp;
        return isl::stat::ok();
      };
      qcopy.foreach_piece(analyzePoly);
      qcopy.release();
      isl::pw_qpolynomial polyModified = clean_poly_lcm(tqp);
      tqp.release();
      auto analyzePolyAgain = [&](isl::set s, isl::qpolynomial qp) {
        tempqp = qp;
        return isl::stat::ok();
      };
      polyModified.foreach_piece(analyzePolyAgain);

      int dims = isl_set_n_dim(trial_set_new);
      isl::set domain = isl::manage(isl_set_copy(trial_set_new));
      std::vector<isl::basic_set> domains;
      auto getDomains = [&](isl::basic_set bset) {
        domains.push_back(bset);
        return isl::stat::ok();
      };
      domain.foreach_basic_set(getDomains);

      for (auto &Domain_c : domains) {

        ver_num = 0;
        isl::basic_set Domain_ = Domain_c.remove_divs();
        if (!interval_bounds) {
          isl_vertices *ver = isl_basic_set_compute_vertices(Domain_.get());
          struct check_poly data3 = {_polynomial, Domain_c, cache_limit};
          isl_vertices_foreach_vertex(ver, &vertex_test, &data3);
          isl_vertices_free(ver);
          if (std::adjacent_find(count_poly.begin(), count_poly.end(),
                                 std::not_equal_to<>()) == count_poly.end()) {
            if (count_poly[0] == false) {
              totalMissCount = 0;
              Domain_.release();
              continue;
            } else {
              isl::val res =
                  isl::manage(isl_set_card(isl::set(Domain_c).release())).max();
              totalMissCount = (long long)(res.get_num_si());
              cache_misses += totalMissCount;
              res.release();
              continue;
            }
          }
        }

        std::vector<isl::qpolynomial> mono;
        auto getMonomial = [&](isl::constraint cons) {
          isl::aff taf = cons.get_aff();
          isl::qpolynomial q =
              isl::manage(isl_qpolynomial_from_aff(taf.release()));
          mono.push_back(q);
          return isl::stat::ok();
        };

        Domain_.foreach_constraint(getMonomial);
        std::vector<isl::qpolynomial> products = BullsEye::_combinations_(
            polyModified.release(), mono.size(), 2, mono);
        products.push_back(tempqp);
        tempqp.release();
        // Get number of rows in matrix
        int combs = BullsEye::_combineNCR_(dims, 2);
        combs += 1 + 2 * dims;
        size_t products_size = products.size();
        // Create Handelman Matrix
        std::vector<std::vector<long>> hmatrix(
            combs, std::vector<long>(products_size, 0));
        int col = 0;
        std::string bmask = "";
        for (auto k : products) {
          isl::qpolynomial k_copy(k);
          auto analyzeTerms = [&](isl::term Term) {
            std::vector<int> Variables(Term.dim(isl::dim::set), 0);
            isl::qpolynomial q;
            for (int i = 0; i < Term.dim(isl::dim::set); i++) {
              Variables[i] += Term.get_exp(isl::dim::set, i);
              int val = Term.get_exp(isl::dim::set, i);
              bmask += std::to_string(val);
              if (val > 1) {
                q = isl::manage(
                    isl_qpolynomial_from_term(isl_term_copy(Term.get())));
                k_copy = isl::manage(
                    isl_qpolynomial_sub(k_copy.release(), q.release()));
              }
            }
            isl::val coeff = Term.get_coefficient_val();
            setcolumn(bmask, col, coeff.get_num_si(), hmatrix);
            coeff.release();
            bmask = "";
            return isl::stat::ok();
          };
          k.foreach_term(analyzeTerms);
          k.release();
          k_copy.release();
          col += 1;
        }

        std::vector<std::vector<long>> transpose(products_size,
                                                 std::vector<long>(combs, 0));

        // Transpose
        for (int m = 0; m < hmatrix.size(); ++m)
          for (int j = 0; j < hmatrix[m].size(); ++j) {
            transpose[j][m] = hmatrix[m][j];
          }

        isl::basic_set working_set =
            isl::manage(isl_basic_set_copy(Domain_c.release()));
        try {
          int limit = cache_size[0];
          isl::space sp = isl::manage(isl_set_get_space(trial_set));
          isl::basic_set init_bset = isl::basic_set::universe(sp);
          isl::local_space ls =
              isl::manage(isl_local_space_from_space(sp.release()));
          // Find Octagon Bounds
          if (!interval_bounds && (isl_set_n_dim(trial_set) == 2)) {
            struct data_poly data = {0,         1,     dims,    products_size,
                                     transpose, limit, hmatrix, ls};
            isl_vertices *verts =
                isl_basic_set_compute_vertices(Domain_.release());
            isl_vertices_foreach_vertex(verts, vertex_count, &data);
            isl_vertices_free(verts);
            for (auto con : constraint_list) {
              init_bset = isl::manage(isl_basic_set_add_constraint(
                  init_bset.release(), con.release()));
            }
          }
          if (interval_bounds && (isl_set_n_dim(trial_set) == 2)) {
            struct data_poly data = {0,         1,     dims,    products_size,
                                     transpose, limit, hmatrix, ls};
            isl_vertices *verts = isl_basic_set_compute_vertices(Domain_.get());
            isl_vertices_foreach_vertex(verts, vertex_count_interval, &data);
            isl_vertices_free(verts);
            isl_vertices *verts2 =
                isl_basic_set_compute_vertices(Domain_.release());
            struct data_poly data2 = {1,         0,     dims,    products_size,
                                      transpose, limit, hmatrix, ls};
            isl_vertices_foreach_vertex(verts2, vertex_count_interval, &data2);
            isl_vertices_free(verts2);
            for (auto con : constraint_list) {
              init_bset = isl::manage(isl_basic_set_add_constraint(
                  init_bset.release(), con.release()));
            }
          }
          if (!interval_bounds && (isl_set_n_dim(trial_set) == 3)) {
            struct data_poly data = {0,         1,     dims,    products_size,
                                     transpose, limit, hmatrix, ls};
            isl_vertices *verts = isl_basic_set_compute_vertices(Domain_.get());
            isl_vertices_foreach_vertex(verts, vertex_count, &data);
            isl_vertices_free(verts);
            isl_vertices *verts2 =
                isl_basic_set_compute_vertices(Domain_.get());
            struct data_poly data2 = {1,         2,     dims,    products_size,
                                      transpose, limit, hmatrix, ls};
            isl_vertices_foreach_vertex(verts2, vertex_count, &data2);
            isl_vertices_free(verts2);
            isl_vertices *verts3 =
                isl_basic_set_compute_vertices(Domain_.release());
            struct data_poly data3 = {2,         0,     dims,    products_size,
                                      transpose, limit, hmatrix, ls};
            isl_vertices_foreach_vertex(verts3, vertex_count, &data3);
            isl_vertices_free(verts3);
            for (auto con : constraint_list) {
              init_bset = isl::manage(isl_basic_set_add_constraint(
                  init_bset.release(), con.release()));
            }
          }
          if (interval_bounds && (isl_set_n_dim(trial_set) == 3)) {
            struct data_poly data = {0,         1,     dims,    products_size,
                                     transpose, limit, hmatrix, ls};
            isl_vertices *verts = isl_basic_set_compute_vertices(Domain_.get());
            isl_vertices_foreach_vertex(verts, vertex_count_interval, &data);
            isl_vertices_free(verts);
            isl_vertices *verts2 =
                isl_basic_set_compute_vertices(Domain_.get());
            struct data_poly data2 = {1,         2,     dims,    products_size,
                                      transpose, limit, hmatrix, ls};
            isl_vertices_foreach_vertex(verts2, vertex_count_interval, &data2);
            isl_vertices_free(verts2);
            isl_vertices *verts3 =
                isl_basic_set_compute_vertices(Domain_.release());
            struct data_poly data3 = {2,         0,     dims,    products_size,
                                      transpose, limit, hmatrix, ls};
            isl_vertices_foreach_vertex(verts3, vertex_count_interval, &data3);
            isl_vertices_free(verts3);
            for (auto con : constraint_list) {
              init_bset = isl::manage(isl_basic_set_add_constraint(
                  init_bset.release(), con.release()));
            }
          }

          working_set = isl::manage(isl_basic_set_intersect(
              working_set.release(), init_bset.release()));
          isl::val res =
              isl::manage(isl_set_card(isl::set(working_set).release())).max();
          totalMissCount = (long long)(res.get_num_si());
          cache_misses += totalMissCount;
          res.release();
          sp.release();
          ls.release();

        } catch (IloException &e) {
          std::cerr << "Error in counting approximate cache misses: " << e
                    << std::endl;
          e.end();
        }

        products.clear();
      }

      Misses.push_back(cache_misses);
      cache_number++;
      domain.release();
    }
  }

  isl_set_free(trial_set);
  isl_set_free(trial_set_new);

  return Misses;
}

isl::aff BullsEye::extractAffineExpression(isl::qpolynomial Polynomial,
                                           isl::set Domain,
                                           std::map<int, long> Values) {
  // buffer the fixed divisors
  std::map<int, isl::aff> Divisors;
  // accumulate the values of all terms
  isl::local_space DS = isl::local_space(Domain.get_space());
  isl::val Zero = isl::val::zero(Domain.get_ctx());
  isl::aff Sum(DS, Zero);
  auto analyzeTerms = [&](isl::term Term) {
    auto Coefficient = Term.get_coefficient_val();
    isl::aff Product(DS, Coefficient);
    // get the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); ++i) {
      int Exponent = Term.get_exp(isl::dim::param, i);
      assert(Exponent <= 1);
      if (Exponent == 1) {
        auto Parameter = isl::aff::var_on_domain(DS, isl::dim::param, i);
        Product = Product.mul(Parameter);
      }
    }
    // get the input exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      int Exponent = Term.get_exp(isl::dim::set, i);
      if (Exponent > 0) {
        if (Values.count(i)) {
          long Coefficient = compute_power(Values[i], Exponent);
          isl::aff Value(DS, isl::val(Domain.get_ctx(), Coefficient));
          Product = Product.mul(Value);
        } else {
          assert(Exponent == 1);
          auto Variable = isl::aff::var_on_domain(DS, isl::dim::set, i);
          Product = Product.mul(Variable);
        }
      }
    }
    // get the divisors exponents
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      int Exponent = Term.get_exp(isl::dim::div, i);
      if (Exponent > 0) {
        // update the divisor if we don't have it yet
        if (Divisors.count(i) == 0) {
          auto Divisor = Term.get_div(i);
          auto Replacement = isl::aff(DS, Divisor.get_constant_val());
          for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
            isl::val Coefficient =
                Divisor.get_coefficient_val(isl::dim::param, j);
            if (!Coefficient.is_zero()) {
              Replacement = Replacement.add_coefficient_val(isl::dim::param, j,
                                                            Coefficient);
            }
          }
          for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
            isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
            if (!Coefficient.is_zero()) {
              if (Values.count(j)) {
                Replacement =
                    Replacement.add_constant_val(Coefficient.mul_ui(Values[j]));
              } else {
                Replacement = Replacement.add_coefficient_val(isl::dim::in, j,
                                                              Coefficient);
              }
            }
          }
          Divisors[i] = Replacement;
        }
        // distinguish cases
        if (Divisors[i].is_cst()) {
          long Base = isl::get_value(Divisors[i].floor().get_constant_val());
          long Coefficient = compute_power(Base, Exponent);
          isl::aff Value(DS, isl::val(Domain.get_ctx(), Coefficient));
          Product = Product.mul(Value);
        } else {
          assert(Exponent == 1);
          Product = Product.mul(Divisors[i].floor());
        }
      }
    }
    Sum = Sum.add(Product);
    return isl::stat::ok();
  };
  Polynomial.foreach_term(analyzeTerms);
  return Sum;
}

isl::pw_aff BullsEye::extractAffineExpression(piece Piece) {
  isl::local_space DS = isl::local_space(Piece.Domain.get_space());
  // initialize the expression with the neutral element
  isl::val Zero = isl::val::zero(Piece.Domain.get_ctx());
  isl::pw_aff Sum(Piece.Domain, Zero);
  for (auto &Term : Piece.Terms) {
    // compute the product
    isl::pw_aff Product(Piece.Domain, Term.Coefficient);
    // multiply with parameter
    for (int i = 0; i < Term.Parameters.size(); ++i) {
      if (Term.Parameters[i] > 0) {
        assert(Term.Parameters[i] <= 1);
        isl::pw_aff Parameter =
            isl::pw_aff::var_on_domain(DS, isl::dim::param, i);
        Product = Product.mul(Parameter);
      }
    }
    // multiply with variable
    for (int i = 0; i < Term.Variables.size(); ++i) {
      if (Term.Variables[i] > 0) {
        assert(Term.Variables[i] <= 1);
        isl::pw_aff Variable;
        Variable = isl::pw_aff::var_on_domain(DS, isl::dim::set, i);
        Product = Product.mul(Variable);
      }
    }
    // multiply with divisor
    auto multiplyTerms = [&](isl::term Term) {
      for (int i = 0; i < Term.dim(isl::dim::div); i++) {
        if (Term.get_exp(isl::dim::div, i) > 0) {
          assert(Term.get_exp(isl::dim::div, i) == 1);
          Product = Product.mul(Term.get_div(i).floor());
        }
      }
      return isl::stat::ok();
    };
    Term.Polynomial.foreach_term(multiplyTerms);
    Sum = Sum.add(Product);
  }
  return Sum;
}

// Naive implementation of setting the coefficients in the matrix
void BullsEye::setcolumn(std::string bmask, int column, long long coeff,
                         std::vector<std::vector<long>> &hm) {

  int row = 0;
  // For Size 2
  if (bmask.size() == 2) {

    if (bmask.compare("00") == 0) {
      row = 0;
      hm[row][column] = coeff;
    } else if (bmask.compare("10") == 0) {
      row = 1;
      hm[row][column] = coeff;
    } else if (bmask.compare("01") == 0) {
      row = 2;
      hm[row][column] = coeff;
    } else if (bmask.compare("11") == 0) {
      row = 3;
      hm[row][column] = coeff;
    } else if (bmask.compare("20") == 0) {
      row = 4;
      hm[row][column] = coeff;
    } else if (bmask.compare("02") == 0) {
      row = 5;
      hm[row][column] = coeff;
    }
    // For Size 1
  } else if (bmask.size() == 1) {
    if (bmask.compare("0") == 0) {
      row = 0;
      hm[row][column] = coeff;
    } else if (bmask.compare("1") == 0) {
      row = 1;
      hm[row][column] = coeff;
    } else if (bmask.compare("2") == 0) {
      row = 2;
      hm[row][column] = coeff;
    }
    // For Size 3
  } else if (bmask.size() == 3) {

    if (bmask.compare("000") == 0) {
      row = 0;
      hm[row][column] = coeff;
    } else if (bmask.compare("100") == 0) {
      row = 1;
      hm[row][column] = coeff;
    } else if (bmask.compare("010") == 0) {
      row = 2;
      hm[row][column] = coeff;
    } else if (bmask.compare("001") == 0) {
      row = 3;
      hm[row][column] = coeff;
    } else if (bmask.compare("110") == 0) {
      row = 4;
      hm[row][column] = coeff;
    } else if (bmask.compare("101") == 0) {
      row = 5;
      hm[row][column] = coeff;
    } else if (bmask.compare("011") == 0) {
      row = 6;
      hm[row][column] = coeff;
    } else if (bmask.compare("200") == 0) {
      row = 7;
      hm[row][column] = coeff;
    } else if (bmask.compare("020") == 0) {
      row = 8;
      hm[row][column] = coeff;
    } else if (bmask.compare("002") == 0) {
      row = 9;
      hm[row][column] = coeff;
    }
    // For Size 6
  } else if (bmask.size() == 6) {

    if (bmask.compare("000000") == 0) {
      row = 0;
      hm[row][column] = coeff;
    } else if (bmask.compare("100000") == 0) {
      row = 1;
      hm[row][column] = coeff;
    } else if (bmask.compare("010000") == 0) {
      row = 2;
      hm[row][column] = coeff;
    } else if (bmask.compare("001000") == 0) {
      row = 3;
      hm[row][column] = coeff;
    } else if (bmask.compare("000100") == 0) {
      row = 4;
      hm[row][column] = coeff;
    } else if (bmask.compare("000010") == 0) {
      row = 5;
      hm[row][column] = coeff;
    } else if (bmask.compare("000001") == 0) {
      row = 6;
      hm[row][column] = coeff;
    } else if (bmask.compare("110000") == 0) {
      row = 7;
      hm[row][column] = coeff;
    } else if (bmask.compare("101000") == 0) {
      row = 8;
      hm[row][column] = coeff;
    } else if (bmask.compare("100100") == 0) {
      row = 9;
      hm[row][column] = coeff;
    } else if (bmask.compare("100010") == 0) {
      row = 10;
      hm[row][column] = coeff;
    } else if (bmask.compare("100001") == 0) {
      row = 11;
      hm[row][column] = coeff;
    } else if (bmask.compare("011000") == 0) {
      row = 12;
      hm[row][column] = coeff;
    } else if (bmask.compare("010100") == 0) {
      row = 13;
      hm[row][column] = coeff;
    } else if (bmask.compare("0100100") == 0) {
      row = 14;
      hm[row][column] = coeff;
    } else if (bmask.compare("010001") == 0) {
      row = 15;
      hm[row][column] = coeff;
    } else if (bmask.compare("001100") == 0) {
      row = 16;
      hm[row][column] = coeff;
    } else if (bmask.compare("001010") == 0) {
      row = 17;
      hm[row][column] = coeff;
    } else if (bmask.compare("001001") == 0) {
      row = 18;
      hm[row][column] = coeff;
    } else if (bmask.compare("000110") == 0) {
      row = 19;
      hm[row][column] = coeff;
    } else if (bmask.compare("000101") == 0) {
      row = 20;
      hm[row][column] = coeff;
    } else if (bmask.compare("000011") == 0) {
      row = 21;
      hm[row][column] = coeff;
    } else if (bmask.compare("200000") == 0) {
      row = 22;
      hm[row][column] = coeff;
    } else if (bmask.compare("020000") == 0) {
      row = 23;
      hm[row][column] = coeff;
    } else if (bmask.compare("002000") == 0) {
      row = 24;
      hm[row][column] = coeff;
    } else if (bmask.compare("000200") == 0) {
      row = 25;
      hm[row][column] = coeff;
    } else if (bmask.compare("000020") == 0) {
      row = 26;
      hm[row][column] = coeff;
    } else if (bmask.compare("000002") == 0) {
      row = 27;
      hm[row][column] = coeff;
    } else if (bmask.compare("111000") == 0) {
      row = 28;
      hm[row][column] = coeff;
    } else if (bmask.compare("110100") == 0) {
      row = 29;
      hm[row][column] = coeff;
    } else if (bmask.compare("110010") == 0) {
      row = 30;
      hm[row][column] = coeff;
    } else if (bmask.compare("110001") == 0) {
      row = 31;
      hm[row][column] = coeff;
    } else if (bmask.compare("101100") == 0) {
      row = 32;
      hm[row][column] = coeff;
    } else if (bmask.compare("101010") == 0) {
      row = 33;
      hm[row][column] = coeff;
    } else if (bmask.compare("101001") == 0) {
      row = 34;
      hm[row][column] = coeff;
    } else if (bmask.compare("100110") == 0) {
      row = 35;
      hm[row][column] = coeff;
    } else if (bmask.compare("100101") == 0) {
      row = 36;
      hm[row][column] = coeff;
    } else if (bmask.compare("100011") == 0) {
      row = 37;
      hm[row][column] = coeff;
    } else if (bmask.compare("210000") == 0) {
      row = 38;
      hm[row][column] = coeff;
    } else if (bmask.compare("201000") == 0) {
      row = 39;
      hm[row][column] = coeff;
    } else if (bmask.compare("200100") == 0) {
      row = 40;
      hm[row][column] = coeff;
    } else if (bmask.compare("200010") == 0) {
      row = 41;
      hm[row][column] = coeff;
    } else if (bmask.compare("200001") == 0) {
      row = 42;
      hm[row][column] = coeff;
    } else if (bmask.compare("120000") == 0) {
      row = 43;
      hm[row][column] = coeff;
    } else if (bmask.compare("021000") == 0) {
      row = 44;
      hm[row][column] = coeff;
    } else if (bmask.compare("020100") == 0) {
      row = 45;
      hm[row][column] = coeff;
    } else if (bmask.compare("020010") == 0) {
      row = 46;
      hm[row][column] = coeff;
    } else if (bmask.compare("020001") == 0) {
      row = 47;
      hm[row][column] = coeff;
    } else if (bmask.compare("102000") == 0) {
      row = 48;
      hm[row][column] = coeff;
    } else if (bmask.compare("012000") == 0) {
      row = 49;
      hm[row][column] = coeff;
    } else if (bmask.compare("002100") == 0) {
      row = 50;
      hm[row][column] = coeff;
    } else if (bmask.compare("002010") == 0) {
      row = 51;
      hm[row][column] = coeff;
    } else if (bmask.compare("002001") == 0) {
      row = 52;
      hm[row][column] = coeff;
    } else if (bmask.compare("100200") == 0) {
      row = 53;
      hm[row][column] = coeff;
    } else if (bmask.compare("010200") == 0) {
      row = 54;
      hm[row][column] = coeff;
    } else if (bmask.compare("001200") == 0) {
      row = 55;
      hm[row][column] = coeff;
    } else if (bmask.compare("000210") == 0) {
      row = 56;
      hm[row][column] = coeff;
    } else if (bmask.compare("000201") == 0) {
      row = 57;
      hm[row][column] = coeff;
    } else if (bmask.compare("100020") == 0) {
      row = 58;
      hm[row][column] = coeff;
    } else if (bmask.compare("010020") == 0) {
      row = 59;
      hm[row][column] = coeff;
    } else if (bmask.compare("001020") == 0) {
      row = 60;
      hm[row][column] = coeff;
    } else if (bmask.compare("000120") == 0) {
      row = 61;
      hm[row][column] = coeff;
    } else if (bmask.compare("000021") == 0) {
      row = 62;
      hm[row][column] = coeff;
    } else if (bmask.compare("100002") == 0) {
      row = 63;
      hm[row][column] = coeff;
    } else if (bmask.compare("010002") == 0) {
      row = 64;
      hm[row][column] = coeff;
    } else if (bmask.compare("001002") == 0) {
      row = 65;
      hm[row][column] = coeff;
    } else if (bmask.compare("000102") == 0) {
      row = 66;
      hm[row][column] = coeff;
    } else if (bmask.compare("000012") == 0) {
      row = 67;
      hm[row][column] = coeff;
    } else if (bmask.compare("300000") == 0) {
      row = 68;
      hm[row][column] = coeff;
    } else if (bmask.compare("030000") == 0) {
      row = 69;
      hm[row][column] = coeff;
    } else if (bmask.compare("003000") == 0) {
      row = 70;
      hm[row][column] = coeff;
    } else if (bmask.compare("000300") == 0) {
      row = 71;
      hm[row][column] = coeff;
    } else if (bmask.compare("000030") == 0) {
      row = 72;
      hm[row][column] = coeff;
    } else if (bmask.compare("000003") == 0) {
      row = 73;
      hm[row][column] = coeff;
    }
  }
}

std::vector<long> BullsEye::countAffineDimensions(piece Piece,
                                                  std::vector<long> Limits) {
  std::vector<long> Results(Limits.size());
  isl::pw_aff LHS = Piece.Expression;
  for (int i = 0; i < Limits.size(); ++i) {
    isl::pw_aff RHS(Piece.Domain, isl::val(Piece.Domain.get_ctx(), Limits[i]));
    isl::set Misses = LHS.gt_set(RHS);
    auto Variable = isl::manage(isl_set_card(Misses.release()));
    assert(isl::get_value(Variable.min()) == isl::get_value(Variable.max()));
    Results[i] = isl::get_value(Variable.max());
  }
  return Results;
}