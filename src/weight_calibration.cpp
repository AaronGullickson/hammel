// Weight-calibration solver declared in weight_calibration.h, plus internal
// inspection entry points used by the tests.
#include "weight_calibration.h"

#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "matching.h"
#include "rate_bands.h"

using namespace Rcpp;

namespace hammel {

// Solve a small dense system with partial pivoting; false if singular.
bool solve_linear(
  std::vector<double> a,
  std::vector<double> b,
  int n,
  std::vector<double>& x
) {
  for (int c = 0; c < n; ++c) {
    int pivot = c;
    for (int r = c + 1; r < n; ++r) {
      if (std::abs(a[r * n + c]) > std::abs(a[pivot * n + c])) {
        pivot = r;
      }
    }
    if (!(std::abs(a[pivot * n + c]) > 1e-300)) {
      return false;
    }
    if (pivot != c) {
      for (int k = 0; k < n; ++k) {
        std::swap(a[c * n + k], a[pivot * n + k]);
      }
      std::swap(b[c], b[pivot]);
    }
    for (int r = c + 1; r < n; ++r) {
      const double factor = a[r * n + c] / a[c * n + c];
      for (int k = c; k < n; ++k) {
        a[r * n + k] -= factor * a[c * n + k];
      }
      b[r] -= factor * b[c];
    }
  }
  x.assign(n, 0);
  for (int r = n - 1; r >= 0; --r) {
    double sum = b[r];
    for (int k = r + 1; k < n; ++k) {
      sum -= a[r * n + k] * x[k];
    }
    x[r] = sum / a[r * n + r];
  }
  return true;
}

WeightSolution WeightCalibrator::solve(
  double now,
  int channel,
  int table_id,
  const std::vector<double>& delta,
  std::vector<double>& warm
) {
  advance(now);
  const LagTable& table = *tables_[table_id];
  const std::vector<double>& sums = sums_[table_id];
  const int G = groups_;
  WeightSolution s;
  s.log_weights.assign(G, 0.0);
  s.target.assign(G, 0.0);
  s.flow.assign(G, 0.0);
  std::vector<double> women_count(G, 0), men_count(G, 0), intensity(G, 0);
  std::vector<std::vector<std::pair<int, double>>> men(G);
  struct Row {
    int group, bin;
    double lambda;
  };
  std::vector<Row> rows;
  for (int g = 0; g < G; ++g) {
    for (int k = 0; k < bins_; ++k) {
      const double age = now - centre(k);
      if (age < cutoff_) {
        continue;
      }
      const size_t i = static_cast<size_t>(g) * bins_ + k;
      if (men_[i] > 0) {
        men[g].push_back({k, men_[i]});
        men_count[g] += men_[i];
      }
      if (women_[i] > 0) {
        women_count[g] += women_[i];
        const double lambda = std::max(0.0, intensities_[i][channel]);
        if (lambda > 0) {
          rows.push_back({g, k, lambda});
          intensity[g] += lambda;
        }
      }
    }
  }
  double total_women = 0, total_intensity = 0;
  for (int g = 0; g < G; ++g) {
    total_women += women_count[g];
    total_intensity += intensity[g];
  }
  if (rows.empty() || !(total_women > 0) || !(total_intensity > 0)) {
    warm = s.log_weights;
    return s;
  }
  // Mirror: men of group g get the per-capita intensity of group-g women
  // (population-wide rate if none are present), then apply target deltas.
  const double overall = total_intensity / total_women;
  std::vector<double> mirror(G, 0);
  double mirror_sum = 0;
  for (int g = 0; g < G; ++g) {
    mirror[g] = men_count[g] *
                (women_count[g] > 0 ? intensity[g] / women_count[g] : overall);
    mirror_sum += mirror[g];
  }
  if (!(mirror_sum > 0)) {
    warm = s.log_weights;
    return s;
  }
  // Normalize log mirrored flows plus delta, preserving zero baseline
  // targets.
  double largest = -std::numeric_limits<double>::infinity();
  for (int g = 0; g < G; ++g) {
    if (mirror[g] > 0) {
      largest = std::max(largest, std::log(mirror[g]) + delta[g]);
    }
  }
  double adjusted_sum = 0;
  for (int g = 0; g < G; ++g) {
    s.target[g] =
      mirror[g] > 0 ? std::exp(std::log(mirror[g]) + delta[g] - largest) : 0;
    adjusted_sum += s.target[g];
  }
  for (int g = 0; g < G; ++g) {
    s.target[g] /= adjusted_sum;
  }
  // Pool sums S[row, male group] come from the maintained sums. Entries below
  // the floor are recomputed exactly; a row that underflows is rescaled in
  // log space.
  const int R = static_cast<int>(rows.size());
  std::vector<double> S(static_cast<size_t>(R) * G, 0);
  // S[r, m], and the lag-table entry of row r against men of group m in bin.
  auto pool = [&](int r, int m) -> double& {
    return S[static_cast<size_t>(r) * G + m];
  };
  auto lag = [&](int r, int m, int bin) {
    return table.index(rows[r].group, m, rows[r].bin - bin + bins_ - 1);
  };
  for (int r = 0; r < R; ++r) {
    double total = 0;
    for (int m = 0; m < G; ++m) {
      // Prohibited pairs contribute exactly zero.
      if (!table.allowed[static_cast<size_t>(rows[r].group) * G + m]) {
        continue;
      }
      double sum = sums
        [(static_cast<size_t>(rows[r].group) * G + m) * bins_ + rows[r].bin];
      if (!(sum >= kPoolSumFloor)) {
        sum = 0;
        for (const auto& bin : men[m]) {
          sum += bin.second * table.weight[lag(r, m, bin.first)];
        }
      }
      pool(r, m) = sum;
      total += sum;
    }
    if (total > 0) {
      continue;
    }
    double maximum = R_NegInf;
    for (int m = 0; m < G; ++m) {
      for (const auto& bin : men[m]) {
        maximum = std::max(maximum, table.utility[lag(r, m, bin.first)]);
      }
    }
    if (maximum == R_NegInf) {
      continue;
    }
    for (int m = 0; m < G; ++m) {
      double sum = 0;
      for (const auto& bin : men[m]) {
        const double u = table.utility[lag(r, m, bin.first)];
        if (u != R_NegInf) {
          sum += bin.second * std::exp(u - maximum);
        }
      }
      pool(r, m) = sum;
    }
  }
  std::vector<double> demand(G, 0);
  for (int r = 0; r < R; ++r) {
    for (int m = 0; m < G; ++m) {
      demand[m] += rows[r].lambda * pool(r, m);
    }
  }
  // Each separate market (D138) is solved on its own: its free groups, with
  // targets normalized within it, over its women's rows. One block is the
  // whole population, as before.
  std::vector<std::vector<int>> block_groups(table.blocks),
    block_rows(table.blocks);
  for (int g = 0; g < G; ++g) {
    if (!(s.target[g] > 0)) {
      s.log_weights[g] = kZeroTargetLogWeight;
    } else if (!(demand[g] > 0)) {
      s.log_weights[g] = kMaxLogWeight;
      s.infeasible = true;
    } else {
      block_groups[table.block[g]].push_back(g);
    }
  }
  for (int r = 0; r < R; ++r) {
    block_rows[table.block[rows[r].group]].push_back(r);
  }
  for (int b = 0; b < table.blocks; ++b) {
    const std::vector<int>& free_groups = block_groups[b];
    const std::vector<int>& block_row_ids = block_rows[b];
    const int n = static_cast<int>(free_groups.size());
    if (n > 1) {
      double target_sum = 0;
      for (int g : free_groups) {
        target_sum += s.target[g];
      }
      std::vector<double> tau(n), x(n, 0);
      for (int j = 0; j < n; ++j) {
        tau[j] = s.target[free_groups[j]] / target_sum;
      }
      if (static_cast<int>(warm.size()) == G) {
        for (int j = 0; j < n; ++j) {
          const double start = warm[free_groups[j]] - warm[free_groups[0]];
          if (R_finite(start) && std::abs(start) < kMaxLogWeight) {
            x[j] = start;
          }
        }
      }
      x[0] = 0;
      std::vector<double> share(n), phi, hess, weights(n);
      // Objective, gradient (phi) and, when `with_hess`, the Hessian. Line
      // search trials need only the objective and gradient.
      auto evaluate = [&](
                        const std::vector<double>& point,
                        std::vector<double>& out_phi,
                        std::vector<double>& out_hess, bool with_hess
                      ) {
        out_phi.assign(n, 0);
        if (with_hess) {
          out_hess.assign(static_cast<size_t>(n) * n, 0);
        }
        double objective = 0, lambda_total = 0;
        for (int j = 0; j < n; ++j) {
          weights[j] = std::exp(point[j]);
        }
        for (int r : block_row_ids) {
          double denominator = 0;
          for (int j = 0; j < n; ++j) {
            share[j] = weights[j] * pool(r, free_groups[j]);
            denominator += share[j];
          }
          if (!(denominator > 0)) {
            continue;
          }
          const double lambda = rows[r].lambda;
          lambda_total += lambda;
          objective += lambda * std::log(denominator);
          for (int j = 0; j < n; ++j) {
            const double p = share[j] / denominator;
            out_phi[j] += lambda * p;
            if (with_hess) {
              for (int k = 0; k < n; ++k) {
                out_hess[static_cast<size_t>(j) * n + k] +=
                  lambda * ((j == k ? p : 0) - p * share[k] / denominator);
              }
            }
          }
        }
        for (double& value : out_phi) {
          value /= lambda_total;
        }
        if (with_hess) {
          for (double& value : out_hess) {
            value /= lambda_total;
          }
        }
        objective /= lambda_total;
        for (int j = 0; j < n; ++j) {
          objective -= tau[j] * point[j];
        }
        return objective;
      };
      bool converged = false;
      std::vector<double> trial(n), trial_phi, trial_hess, step;
      for (int iteration = 0; iteration < 100; ++iteration) {
        const double objective = evaluate(x, phi, hess, true);
        double error = 0;
        for (int j = 0; j < n; ++j) {
          error = std::max(error, std::abs(phi[j] - tau[j]));
        }
        if (error < 1e-10) {
          converged = true;
          break;
        }
        // Newton step over the free log weights (the first is the reference).
        const int m = n - 1;
        std::vector<double> a(static_cast<size_t>(m) * m), b(m);
        for (int j = 0; j < m; ++j) {
          b[j] = -(phi[j + 1] - tau[j + 1]);
          for (int k = 0; k < m; ++k) {
            a[static_cast<size_t>(j) * m + k] =
              hess[static_cast<size_t>(j + 1) * n + k + 1];
          }
          a[static_cast<size_t>(j) * m + j] += 1e-14;
        }
        if (!solve_linear(a, b, m, step)) {
          break;
        }
        // Accept a step that lowers the convex objective or the target
        // mismatch; near the solution, roundoff can hide the objective's
        // decrease.
        double length = 1;
        bool moved = false;
        while (length > 1e-12) {
          trial[0] = 0;
          for (int j = 0; j < m; ++j) {
            trial[j + 1] = std::max(
              -kMaxLogWeight,
              std::min(kMaxLogWeight, x[j + 1] + length * step[j])
            );
          }
          const double trial_objective =
            evaluate(trial, trial_phi, trial_hess, false);
          double trial_error = 0;
          for (int j = 0; j < n; ++j) {
            trial_error =
              std::max(trial_error, std::abs(trial_phi[j] - tau[j]));
          }
          if (trial_objective < objective || trial_error < error) {
            moved = true;
            break;
          }
          length /= 2;
        }
        if (!moved) {
          break;
        }
        x = trial;
      }
      if (!converged) {
        s.infeasible = true;
      }
      for (int j = 0; j < n; ++j) {
        s.log_weights[free_groups[j]] = x[j];
      }
    }
  }
  // Expected flows under the final weights, in log space so that women whose
  // only allowed men carry near-zero weights still count, as they do in
  // actual selection. Women with no allowed men cannot be matched by any
  // weights and are left out. A remaining miss means the target was not met.
  double flow_total = 0;
  std::vector<double> share(G);
  for (int r = 0; r < R; ++r) {
    double top = R_NegInf;
    for (int m = 0; m < G; ++m) {
      const double men_sum = pool(r, m);
      if (men_sum > 0) {
        top = std::max(top, s.log_weights[m] + std::log(men_sum));
      }
    }
    if (top == R_NegInf) {
      continue;
    }
    double denominator = 0;
    for (int m = 0; m < G; ++m) {
      const double men_sum = pool(r, m);
      share[m] =
        men_sum > 0 ? std::exp(s.log_weights[m] + std::log(men_sum) - top) : 0;
      denominator += share[m];
    }
    flow_total += rows[r].lambda;
    for (int m = 0; m < G; ++m) {
      s.flow[m] += rows[r].lambda * share[m] / denominator;
    }
  }
  if (flow_total > 0) {
    for (int m = 0; m < G; ++m) {
      s.flow[m] /= flow_total;
    }
    // With separate markets, each block's targets are rescaled to its share
    // of the flow (its women's searching), which is all it can receive.
    if (table.blocks > 1) {
      std::vector<double> block_flow(table.blocks, 0),
        block_target(table.blocks, 0);
      for (int m = 0; m < G; ++m) {
        block_flow[table.block[m]] += s.flow[m];
        block_target[table.block[m]] += s.target[m];
      }
      for (int m = 0; m < G; ++m) {
        const int b = table.block[m];
        s.target[m] = block_target[b] > 0
                        ? s.target[m] * block_flow[b] / block_target[b]
                        : 0;
      }
    }
    double miss = 0;
    for (int m = 0; m < G; ++m) {
      miss = std::max(miss, std::abs(s.flow[m] - s.target[m]));
    }
    if (miss > 1e-8) {
      s.infeasible = true;
    }
  }
  warm = s.log_weights;
  return s;
}

}  // namespace hammel

using namespace hammel;

// Internal entry point: calibrate weights for a supplied unpartnered state.
// [[Rcpp::export]]
List inspect_weighted_calibration_cpp(
  NumericVector female_born,
  IntegerVector female_groups,
  LogicalVector female_fertile,
  NumericVector male_born,
  IntegerVector male_groups,
  DataFrame rates,
  List matching,
  int group_count,
  double now,
  double cutoff,
  double width,
  bool fertility,
  Nullable<IntegerVector> female_previous = R_NilValue,
  Nullable<NumericVector> female_previous_ends = R_NilValue
) {
  if (female_groups.size() != female_born.size() ||
      female_fertile.size() != female_born.size() ||
      male_groups.size() != male_born.size() || group_count < 1) {
    stop("Invalid weighted-calibration inputs.");
  }
  double lo = now, hi = now;
  for (double born : female_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  for (double born : male_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  const bool previous_specific =
    rates.containsElementNamed("previous_union") && !fertility;
  IntegerVector previous(female_born.size(), 0);
  if (female_previous.isNotNull()) {
    previous = IntegerVector(female_previous);
  }
  if (previous.size() != female_born.size()) {
    stop("Invalid previous-union calibration inputs.");
  }
  for (int p : previous) {
    if (p < 0 || p > 2) {
      stop("Invalid previous-union calibration category.");
    }
  }
  WeightCalibrator calibrator(lo, hi, group_count, width, cutoff);
  const LagTable table =
    build_lag_table(matching, group_count, calibrator.bins(), width);
  calibrator.add_table(&table);
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    calibrator.add_man(male_groups[i], male_born[i], 1);
  }
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    calibrator.add_woman(female_groups[i], female_born[i], 1);
  }
  NumericVector previous_ends(female_born.size(), NA_REAL);
  if (female_previous_ends.isNotNull()) {
    previous_ends = NumericVector(female_previous_ends);
  }
  if (previous_ends.size() != female_born.size()) {
    stop("Invalid previous ending times.");
  }
  const int categories = previous_specific ? 3 : 1;
  std::vector<bool> duration_bands(group_count * categories, false);
  std::vector<std::vector<Band>> bands(group_count * categories);
  std::vector<const std::vector<Band>*> hazards(group_count * categories);
  const std::array<std::string, 3> labels = {"none", "separated", "widowed"};
  for (int g = 0; g < group_count; ++g) {
    for (int p = 0; p < categories; ++p) {
      const int i = g * categories + p;
      duration_bands[i] =
        previous_specific && search_uses_duration(rates, labels[p]);
      bands[i] = read_bands(
        rates, "", "", g, previous_specific ? labels[p] : "", duration_bands[i]
      );
      hazards[i] = &bands[i];
    }
  }
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    const int schedule =
      female_groups[i] * categories + (previous_specific ? previous[i] : 0);
    const double origin =
      duration_bands[schedule] ? previous_ends[i] : female_born[i];
    if (!R_finite(origin)) {
      stop("Duration calibration requires a previous ending time.");
    }
    const double value = current_hazard(hazards[schedule], origin, now).first;
    calibrator.replace_intensity(
      female_groups[i], female_born[i], {0, 0},
      {value, female_fertile[i] ? value : 0}
    );
  }
  std::vector<double> warm;
  const WeightSolution s = calibrator.solve(
    now, fertility, 0, as<std::vector<double>>(matching["delta"]), warm
  );
  return List::create(
    _["log_weights"] = wrap(s.log_weights), _["target"] = wrap(s.target),
    _["flow"] = wrap(s.flow), _["infeasible"] = s.infeasible
  );
}

// Internal entry point: apply a sequence of male pool changes at increasing
// times, then compare incremental pool sums and weights with a full rebuild.
// [[Rcpp::export]]
List inspect_weighted_incremental_cpp(
  NumericVector female_born,
  IntegerVector female_groups,
  NumericVector male_born,
  IntegerVector male_groups,
  NumericVector op_times,
  IntegerVector op_men,
  NumericVector op_deltas,
  DataFrame rates,
  List matching,
  int group_count,
  double cutoff
) {
  if (op_times.size() == 0 || op_men.size() != op_times.size() ||
      op_deltas.size() != op_times.size()) {
    stop("Invalid incremental-calibration inputs.");
  }
  double lo = op_times[0], hi = op_times[op_times.size() - 1];
  for (double born : female_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  for (double born : male_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  WeightCalibrator incremental(lo, hi, group_count, kWeightBinWidth, cutoff);
  const LagTable table =
    build_lag_table(matching, group_count, incremental.bins(), kWeightBinWidth);
  incremental.add_table(&table);
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    incremental.add_woman(female_groups[i], female_born[i], 1);
  }
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    incremental.add_man(male_groups[i], male_born[i], 1);
  }
  std::vector<std::vector<Band>> bands(group_count);
  std::vector<const std::vector<Band>*> hazards(group_count);
  for (int g = 0; g < group_count; ++g) {
    bands[g] = read_bands(rates, "", "", g);
    hazards[g] = &bands[g];
  }
  const std::vector<double> delta = as<std::vector<double>>(matching["delta"]);
  std::vector<double> warm;
  std::vector<double> prior(female_born.size(), 0);
  auto update_intensities = [&](double now) {
    for (R_xlen_t j = 0; j < female_born.size(); ++j) {
      const double value =
        current_hazard(hazards[female_groups[j]], female_born[j], now).first;
      incremental.replace_intensity(
        female_groups[j], female_born[j], {prior[j], prior[j]}, {value, value}
      );
      prior[j] = value;
    }
  };
  for (R_xlen_t i = 0; i < op_times.size(); ++i) {
    if (i % 1024 == 0) {
      checkUserInterrupt();
    }
    update_intensities(op_times[i]);
    incremental.solve(op_times[i], false, 0, delta, warm);
    incremental.add_man(
      male_groups[op_men[i]], male_born[op_men[i]], op_deltas[i]
    );
  }
  const double end = op_times[op_times.size() - 1];
  WeightCalibrator rebuilt = incremental;
  rebuilt.refresh();
  const double difference = incremental.max_sum_difference(rebuilt, 0);
  std::vector<double> warm_a, warm_b;
  const WeightSolution a = incremental.solve(end, false, 0, delta, warm_a);
  const WeightSolution b = rebuilt.solve(end, false, 0, delta, warm_b);
  return List::create(
    _["max_sum_difference"] = difference,
    _["incremental"] = wrap(a.log_weights), _["rebuilt"] = wrap(b.log_weights)
  );
}
