// Matching utilities and exact partner draws declared in matching.h, plus
// internal inspection entry points used by the tests.
#include "matching.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

#include "event_queues.h"

using namespace Rcpp;

namespace hammel {

NumericVector matching_utilities(
  const NumericVector& gaps,
  const List& matching,
  int female_group,
  IntegerVector male_groups,
  IntegerVector female_groups
) {
  NumericVector quadratic = matching["quadratic"];
  NumericVector utility(gaps.size());
  if (quadratic.size() == 2) {
    const bool has_squared = as<bool>(matching["has_squared"]);
    for (R_xlen_t j = 0; j < gaps.size(); ++j) {
      const double squared = has_squared ? gaps[j] * gaps[j] : 0;
      if (!R_finite(squared)) {
        stop(
          "Matching formula must produce finite values for every candidate."
        );
      }
      utility[j] = quadratic[0] * gaps[j] + quadratic[1] * squared;
    }
  } else {
    Function evaluate = matching["evaluate"];
    if (as<bool>(matching["group_formula"])) {
      if ((female_group < 0 && female_groups.size() != gaps.size()) ||
          male_groups.size() != gaps.size()) {
        stop("Matching group predictors require group data.");
      }
      if (female_groups.size()) {
        utility = as<NumericVector>(evaluate(gaps, female_groups, male_groups));
      } else {
        utility = as<NumericVector>(evaluate(gaps, female_group, male_groups));
      }
    } else {
      utility = as<NumericVector>(evaluate(gaps));
    }
  }
  if (utility.size() != gaps.size()) {
    stop("Matching formula returned the wrong number of utilities.");
  }
  NumericMatrix offsets = matching["pair_offsets"];
  if (offsets.nrow()) {
    if ((female_groups.size() != gaps.size() &&
         (female_group < 0 || female_group >= offsets.nrow())) ||
        male_groups.size() != gaps.size()) {
      stop("Pair effects require group data.");
    }
    for (R_xlen_t j = 0; j < utility.size(); ++j) {
      if (male_groups[j] < 0 || male_groups[j] >= offsets.ncol()) {
        stop("Invalid matching group index.");
      }
      const int fg = female_groups.size() ? female_groups[j] : female_group;
      if (fg < 0 || fg >= offsets.nrow()) {
        stop("Invalid matching group index.");
      }
      utility[j] += offsets(fg, male_groups[j]);
    }
  }
  for (R_xlen_t j = 0; j < utility.size(); ++j) {
    if (!R_finite(gaps[j]) || !R_finite(utility[j])) {
      stop("Matching ages and utilities must be finite.");
    }
  }
  return utility;
}

LagTable
build_lag_table(const List& model, int groups, int bins, double width) {
  LagTable table;
  table.groups = groups;
  table.span = 2 * bins - 1;
  const size_t n = static_cast<size_t>(groups) * groups * table.span;
  NumericVector gaps(n);
  IntegerVector female(n), male(n);
  for (int f = 0; f < groups; ++f) {
    for (int m = 0; m < groups; ++m) {
      for (int lag = 0; lag < table.span; ++lag) {
        const size_t i = table.index(f, m, lag);
        gaps[i] = (lag - (bins - 1)) * width;
        female[i] = f;
        male[i] = m;
      }
    }
  }
  const NumericVector utility =
    matching_utilities(gaps, model, -1, male, female);
  const bool prevent_all = as<bool>(model["prevent_exogamy_all"]);
  LogicalMatrix blocked = model["blocked_pairs"];
  table.utility.assign(utility.begin(), utility.end());
  // Excluded pairs get -Inf utility; the allowed pairs connect groups into
  // separate markets (D138), found by union-find.
  table.allowed.assign(static_cast<size_t>(groups) * groups, 0);
  std::vector<int> parent(groups);
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> root = [&](int g) {
    return parent[g] == g ? g : parent[g] = root(parent[g]);
  };
  double maximum = R_NegInf;
  for (int f = 0; f < groups; ++f) {
    for (int m = 0; m < groups; ++m) {
      const bool excluded =
        (prevent_all && f != m) || (blocked.nrow() && blocked(f, m));
      if (!excluded) {
        table.allowed[static_cast<size_t>(f) * groups + m] = 1;
        parent[root(f)] = root(m);
      }
      for (int lag = 0; lag < table.span; ++lag) {
        double& value = table.utility[table.index(f, m, lag)];
        if (excluded) {
          value = R_NegInf;
        } else {
          maximum = std::max(maximum, value);
        }
      }
    }
  }
  table.weight.resize(n);
  for (size_t i = 0; i < n; ++i) {
    table.weight[i] =
      table.utility[i] == R_NegInf ? 0 : std::exp(table.utility[i] - maximum);
  }
  table.block.assign(groups, 0);
  std::vector<int> label(groups, -1);
  table.blocks = 0;
  for (int g = 0; g < groups; ++g) {
    const int r = root(g);
    if (label[r] < 0) {
      label[r] = table.blocks++;
    }
    table.block[g] = label[r];
  }
  return table;
}

DrawModel make_draw_model(const List& model, int groups, double gap_range) {
  DrawModel draw;
  draw.model = model;
  draw.groups = groups;
  draw.offsets = as<NumericMatrix>(model["pair_offsets"]);
  draw.blocked = as<LogicalMatrix>(model["blocked_pairs"]);
  draw.prevent_all = as<bool>(model["prevent_exogamy_all"]);
  NumericVector q = model["quadratic"];
  draw.native = q.size() == 2;
  if (draw.native) {
    draw.linear = q[0];
    draw.quadratic = as<bool>(model["has_squared"]) ? q[1] : 0.0;
    return draw;
  }
  // Fine grid over every possible gap in the run: 1/32-year resolution.
  draw.grid_step = kWeightBinWidth / 8;
  const double reach = gap_range + 2 * kWeightBinWidth;
  draw.grid_start = -reach;
  draw.grid_size = static_cast<int>(std::ceil(2 * reach / draw.grid_step)) + 1;
  const size_t n = static_cast<size_t>(groups) * groups * draw.grid_size;
  NumericVector gaps(n);
  IntegerVector female(n), male(n);
  for (int f = 0; f < groups; ++f) {
    for (int m = 0; m < groups; ++m) {
      for (int j = 0; j < draw.grid_size; ++j) {
        const size_t i =
          (static_cast<size_t>(f) * groups + m) * draw.grid_size + j;
        gaps[i] = draw.grid_start + j * draw.grid_step;
        female[i] = f;
        male[i] = m;
      }
    }
  }
  const NumericVector values =
    matching_utilities(gaps, model, -1, male, female);
  draw.grid.assign(values.begin(), values.end());
  return draw;
}

int draw_partner(
  const MalePool& pool,
  const DrawModel& draw,
  int female_group,
  double female_born,
  double now,
  double cutoff,
  const std::vector<double>& log_weights,
  const std::function<double(int)>& born_of,
  const std::function<bool(int)>& excluded,
  DrawWorkspace& ws,
  double& approximate
) {
  if (pool.empty()) {
    return -1;
  }
  const int last = pool.last_bin(now, cutoff);
  const double w = pool.width();
  ws.group.clear();
  ws.bin.clear();
  ws.top.clear();
  ws.cumulative.clear();
  double peak = R_NegInf;
  for (int g = 0; g < pool.groups(); ++g) {
    if (!draw.allowed(female_group, g)) {
      continue;
    }
    for (int k = 0; k <= last; ++k) {
      const std::vector<int>& members = pool.members(g, k);
      if (members.empty()) {
        continue;
      }
      const double oldest_birth = static_cast<double>(pool.first() + k) * w;
      const double top = draw.bound(
        female_group, g, female_born - oldest_birth - w,
        female_born - oldest_birth
      );
      const double value =
        log_weights[g] + top + std::log(static_cast<double>(members.size()));
      ws.group.push_back(g);
      ws.bin.push_back(k);
      ws.top.push_back(top);
      ws.cumulative.push_back(value);
      peak = std::max(peak, value);
    }
  }
  if (ws.group.empty() || peak == R_NegInf) {
    return -1;
  }
  double total = 0;
  for (double& value : ws.cumulative) {
    total += std::exp(value - peak);
    value = total;
  }
  for (int attempt = 0; attempt < 200; ++attempt) {
    const double u = R::runif(0, 1) * total;
    const size_t c = std::min(
      static_cast<size_t>(
        std::upper_bound(ws.cumulative.begin(), ws.cumulative.end(), u) -
        ws.cumulative.begin()
      ),
      ws.group.size() - 1
    );
    const std::vector<int>& members = pool.members(ws.group[c], ws.bin[c]);
    const int man = members[uniform_index(R::runif(0, 1), members.size())];
    const double utility =
      draw.utility(female_group, ws.group[c], female_born - born_of(man));
    const double ratio = std::exp(utility - ws.top[c]);
    if (ratio > 1 + 1e-9) {
      approximate += 1;  // custom-formula bound was too low
    }
    if (R::runif(0, 1) >= ratio) {
      continue;
    }
    if (now - born_of(man) < cutoff || excluded(man)) {
      continue;
    }
    return man;
  }
  // Many rejections (for example, most nearby men are kin): enumerate the
  // eligible pool and draw exactly from it.
  std::vector<int> candidates;
  std::vector<double> scores;
  for (int g = 0; g < pool.groups(); ++g) {
    if (!draw.allowed(female_group, g)) {
      continue;
    }
    std::vector<int> men;
    for (int k = 0; k <= last; ++k) {
      for (int man : pool.members(g, k)) {
        if (now - born_of(man) >= cutoff && !excluded(man)) {
          men.push_back(man);
        }
      }
    }
    if (men.empty()) {
      continue;
    }
    NumericVector gaps(men.size());
    for (size_t j = 0; j < men.size(); ++j) {
      gaps[j] = female_born - born_of(men[j]);
    }
    const NumericVector utility = draw.utilities(female_group, g, gaps);
    for (size_t j = 0; j < men.size(); ++j) {
      candidates.push_back(men[j]);
      scores.push_back(log_weights[g] + utility[j]);
    }
  }
  if (candidates.empty()) {
    return -1;
  }
  const double best = *std::max_element(scores.begin(), scores.end());
  double sum = 0;
  for (double& score : scores) {
    sum += std::exp(score - best);
    score = sum;
  }
  const double u = R::runif(0, 1) * sum;
  const size_t j = std::min(
    static_cast<size_t>(
      std::upper_bound(scores.begin(), scores.end(), u) - scores.begin()
    ),
    candidates.size() - 1
  );
  return candidates[j];
}

}  // namespace hammel

using namespace hammel;

// Internal analytical diagnostics use the same native/R scoring kernel.
// [[Rcpp::export]]
NumericVector
inspect_matching_cpp(NumericVector age_difference, List matching) {
  return matching_utilities(age_difference, matching);
}

// Group indices refer to the labels supplied to prepare_matching_model().
// [[Rcpp::export]]
NumericVector inspect_group_matching_cpp(
  NumericVector age_difference,
  List matching,
  int female_group,
  IntegerVector male_groups
) {
  return matching_utilities(
    age_difference, matching, female_group, male_groups
  );
}

// Internal entry point: exact draws from a frozen pool, with the exact
// full-pool conditional-logit probabilities for comparison.
// [[Rcpp::export]]
List inspect_exact_draw_cpp(
  NumericVector male_born,
  IntegerVector male_group,
  LogicalVector male_excluded,
  double female_born,
  int female_group,
  NumericVector log_weights,
  List matching,
  int groups,
  double now,
  double cutoff,
  int draws
) {
  if (male_group.size() != male_born.size() ||
      male_excluded.size() != male_born.size() ||
      log_weights.size() != groups || draws < 0) {
    stop("Invalid exact-draw inputs.");
  }
  double lo = std::min(female_born, now), hi = now;
  for (double born : male_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  MalePool pool(lo, hi, groups, kWeightBinWidth);
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    pool.add(i, male_group[i], male_born[i]);
  }
  const DrawModel draw = make_draw_model(matching, groups, hi - lo);
  const std::vector<double> weights(log_weights.begin(), log_weights.end());
  auto born_of = [&](int man) { return male_born[man]; };
  auto excluded = [&](int man) {
    return static_cast<bool>(male_excluded[man]);
  };
  DrawWorkspace ws;
  double approximate = 0;
  IntegerVector counts(male_born.size());
  int failures = 0;
  for (int d = 0; d < draws; ++d) {
    if (d % 4096 == 0) {
      checkUserInterrupt();
    }
    const int man = draw_partner(
      pool, draw, female_group, female_born, now, cutoff, weights, born_of,
      excluded, ws, approximate
    );
    if (man < 0) {
      ++failures;
    } else {
      ++counts[man];
    }
  }
  NumericVector score(male_born.size(), R_NegInf);
  double best = R_NegInf;
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    const int g = male_group[i];
    if (now - male_born[i] < cutoff || male_excluded[i] ||
        !draw.allowed(female_group, g)) {
      continue;
    }
    score[i] =
      weights[g] + draw.utility(female_group, g, female_born - male_born[i]);
    best = std::max(best, score[i]);
  }
  NumericVector probability(male_born.size());
  double total = 0;
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    if (score[i] != R_NegInf) {
      probability[i] = std::exp(score[i] - best);
      total += probability[i];
    }
  }
  if (total > 0) {
    for (double& p : probability) {
      p /= total;
    }
  }
  return List::create(
    _["counts"] = counts, _["probability"] = probability,
    _["failures"] = failures, _["approximate"] = approximate
  );
}
