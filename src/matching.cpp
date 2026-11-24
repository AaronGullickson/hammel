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
          "[matching_utility] Matching formula must produce finite values for "
          "every candidate."
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
      stop("[matching_utility] Matching ages and utilities must be finite.");
    }
  }
  return utility;
}

MatchingTypes::MatchingTypes(
  int categories,
  std::vector<int> offsets,
  std::vector<int> sizes,
  int units
)
    : categories_(categories),
      units_(units),
      offsets_(std::move(offsets)),
      sizes_(std::move(sizes)) {
  if (categories < 1 || offsets_.size() != sizes_.size() || units < 0 ||
      (units > 0) != !offsets_.empty()) {
    stop("Invalid matching types.");
  }
  // Without ancestry endogamy, category c is type c with the empty profile.
  if (!active()) {
    profiles_.emplace_back();
    profile_ids_[profiles_.back()] = 0;
    for (int c = 0; c < categories; ++c) {
      type(c, 0);
    }
  }
}

int MatchingTypes::profile_of(const double* shares) {
  if (!active()) {
    return 0;
  }
  std::vector<int> key;
  for (size_t d = 0; d < offsets_.size(); ++d) {
    const std::vector<int> part =
      round_profile(shares + offsets_[d], sizes_[d], units_);
    key.insert(key.end(), part.begin(), part.end());
  }
  const auto found = profile_ids_.find(key);
  if (found != profile_ids_.end()) {
    return found->second;
  }
  const int id = static_cast<int>(profiles_.size());
  profile_ids_[key] = id;
  profiles_.push_back(std::move(key));
  return id;
}

int MatchingTypes::type(int category, int profile) {
  if (category < 0 || category >= categories_ || profile < 0 ||
      profile >= static_cast<int>(profiles_.size())) {
    stop("Invalid matching type.");
  }
  const auto key = std::make_pair(category, profile);
  const auto found = type_ids_.find(key);
  if (found != type_ids_.end()) {
    return found->second;
  }
  const int id = count();
  type_ids_[key] = id;
  category_.push_back(category);
  profile_.push_back(profile);
  return id;
}

std::vector<int> round_profile(const double* shares, int size, int units) {
  std::vector<int> out(size);
  std::vector<double> remainder(size);
  int left = units;
  for (int j = 0; j < size; ++j) {
    const double scaled = shares[j] * units;
    // A small tolerance keeps shares on the grid (0.375 at 1/8) whole.
    out[j] = static_cast<int>(std::floor(scaled + 1e-9));
    remainder[j] = scaled - out[j];
    left -= out[j];
  }
  // Remainders within 1e-9 of each other tie; the earlier component wins.
  std::vector<char> used(size, 0);
  for (; left > 0; --left) {
    int best = -1;
    for (int j = 0; j < size; ++j) {
      if (!used[j] && (best < 0 || remainder[j] > remainder[best] + 1e-9)) {
        best = j;
      }
    }
    if (best < 0) {
      stop("Ancestry shares must sum to one.");
    }
    used[best] = 1;
    ++out[best];
  }
  if (left < 0) {
    stop("Ancestry shares must sum to one.");
  }
  return out;
}

double
AncestryTerm::value(const MatchingTypes& types, int female, int male) const {
  if (lambda.empty()) {
    return 0;
  }
  const std::vector<int>& a = types.profile_units(female);
  const std::vector<int>& b = types.profile_units(male);
  double sum = 0;
  for (size_t k = 0; k < lambda.size(); ++k) {
    sum += lambda[k] * std::abs(a[k] - b[k]);
  }
  return -sum / types.units();
}

double
AncestryTerm::factor(const MatchingTypes& types, int female, int male) const {
  if (lambda.empty()) {
    return 1;
  }
  if (female >= static_cast<int>(cache.size())) {
    cache.resize(female + 1);
  }
  std::vector<double>& row = cache[female];
  if (male >= static_cast<int>(row.size())) {
    row.resize(male + 1, R_NaN);
  }
  if (std::isnan(row[male])) {
    row[male] = std::exp(value(types, female, male));
  }
  return row[male];
}

AncestryTerm read_ancestry_term(const List& model) {
  AncestryTerm term;
  if (model.containsElementNamed("ancestry_lambda")) {
    const NumericVector lambda = model["ancestry_lambda"];
    term.lambda.assign(lambda.begin(), lambda.end());
  }
  return term;
}

LagTable
build_lag_table(const List& model, int groups, int bins, double width) {
  LagTable table;
  table.groups = groups;
  table.ancestry = read_ancestry_term(model);
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
  draw.ancestry = read_ancestry_term(model);
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
  const MatchingTypes& types,
  int female_type,
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
  const int female_group = types.category(female_type);
  const int female_profile = types.profile(female_type);
  // The ancestry term is constant within a cell, so bounds stay exact.
  auto ancestry = [&](int type) {
    return draw.ancestry.value(types, female_profile, types.profile(type));
  };
  const int last = pool.last_bin(now, cutoff);
  const double w = pool.width();
  ws.group.clear();
  ws.bin.clear();
  ws.top.clear();
  ws.cumulative.clear();
  double peak = R_NegInf;
  for (int g = 0; g < pool.groups(); ++g) {
    const int male_group = types.category(g);
    if (!draw.allowed(female_group, male_group)) {
      continue;
    }
    const double term = ancestry(g);
    for (int k = 0; k <= last; ++k) {
      const std::vector<int>& members = pool.members(g, k);
      if (members.empty()) {
        continue;
      }
      const double oldest_birth = static_cast<double>(pool.first() + k) * w;
      const double top =
        draw.bound(
          female_group, male_group, female_born - oldest_birth - w,
          female_born - oldest_birth
        ) +
        term;
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
    const int type = ws.group[c];
    const double utility =
      draw.utility(
        female_group, types.category(type), female_born - born_of(man)
      ) +
      ancestry(type);
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
    const int male_group = types.category(g);
    if (!draw.allowed(female_group, male_group)) {
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
    const NumericVector utility =
      draw.utilities(female_group, male_group, gaps);
    const double term = ancestry(g);
    for (size_t j = 0; j < men.size(); ++j) {
      candidates.push_back(men[j]);
      scores.push_back(log_weights[g] + utility[j] + term);
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
// full-pool conditional-logit probabilities for comparison. With an ancestry
// layout, rows of `male_shares` and `female_shares` give ancestry shares and
// `log_weights` are by category (shared by a category's types).
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
  int draws,
  Nullable<List> ancestry_layout = R_NilValue,
  Nullable<NumericMatrix> male_shares = R_NilValue,
  Nullable<NumericVector> female_shares = R_NilValue
) {
  if (male_group.size() != male_born.size() ||
      male_excluded.size() != male_born.size() ||
      log_weights.size() != groups || draws < 0) {
    stop("Invalid exact-draw inputs.");
  }
  MatchingTypes types = inspect_types(groups, ancestry_layout);
  std::vector<int> male_type(male_born.size());
  int female_type = female_group;
  if (types.active()) {
    const NumericMatrix shares(male_shares.get());
    const NumericVector own(female_shares.get());
    for (R_xlen_t i = 0; i < male_born.size(); ++i) {
      const NumericVector row = shares(i, _);
      male_type[i] = types.type(male_group[i], types.profile_of(row.begin()));
    }
    female_type = types.type(female_group, types.profile_of(own.begin()));
  } else {
    std::copy(male_group.begin(), male_group.end(), male_type.begin());
  }
  double lo = std::min(female_born, now), hi = now;
  for (double born : male_born) {
    lo = std::min(lo, born);
    hi = std::max(hi, born);
  }
  MalePool pool(lo, hi, groups, kWeightBinWidth);
  pool.grow(types.count());
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    pool.add(i, male_type[i], male_born[i]);
  }
  const DrawModel draw = make_draw_model(matching, groups, hi - lo);
  std::vector<double> weights(types.count());
  for (int t = 0; t < types.count(); ++t) {
    weights[t] = log_weights[types.category(t)];
  }
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
      pool, draw, types, female_type, female_born, now, cutoff, weights,
      born_of, excluded, ws, approximate
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
    score[i] = weights[male_type[i]] +
               draw.utility(female_group, g, female_born - male_born[i]) +
               draw.ancestry.value(
                 types, types.profile(female_type), types.profile(male_type[i])
               );
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

// Internal entry point: each row's shares rounded to whole units (D199).
// [[Rcpp::export]]
IntegerMatrix inspect_round_profile_cpp(NumericMatrix shares, int units) {
  IntegerMatrix out(shares.nrow(), shares.ncol());
  for (int i = 0; i < shares.nrow(); ++i) {
    const NumericVector row = shares(i, _);
    const std::vector<int> rounded =
      round_profile(row.begin(), shares.ncol(), units);
    for (int j = 0; j < shares.ncol(); ++j) {
      out(i, j) = rounded[j];
    }
  }
  return out;
}
