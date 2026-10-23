#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

using namespace Rcpp;

// One piecewise constant age interval: [start, end), in years, with an annual
// person-year rate. An open final interval has end = infinity.
struct Band {
  double start, end, rate;
};

// Copy the R rate table into a growable C++ collection. Mortality selects
// one sex; fertility uses all rows. Sort the copy without changing R inputs.
std::vector<Band> read_bands(DataFrame rates, const std::string& sex = "",
                             const std::string& status = "", int group = -1) {
  NumericVector starts = rates["age_start"], ends = rates["age_end"];
  NumericVector values = rates["rate"];
  CharacterVector sexes;
  if (!sex.empty()) sexes = rates["sex"];
  CharacterVector statuses;
  if (!status.empty()) statuses = rates["union_status"];
  IntegerVector groups;
  const bool grouped = rates.containsElementNamed(".hammel_group");
  if (grouped) groups = rates[".hammel_group"];
  std::vector<Band> bands;
  for (int i = 0; i < starts.size(); ++i) {
    if ((sex.empty() || as<std::string>(sexes[i]) == sex) &&
        (status.empty() || as<std::string>(statuses[i]) == status) &&
        (!grouped || groups[i] == group))
      bands.push_back({starts[i], ends[i], values[i]});
  }
  std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) {
    return a.start < b.start;
  });
  return bands;
}

// Calculate waiting time until the next event. The threshold is drawn
// from Exp(1), equivalent in distribution to -log(U), with U uniform
// on (0, 1). Subtract rate * duration across successive age bands.
// In the band where the remaining threshold is reached, divide it
// by the rate to obtain the remaining waiting time.
double waiting_time(const std::vector<Band>& bands, double age, double threshold) {
  const double origin = age;
  for (const Band& band : bands) {
    // if we are older than the end of this band, skip it
    if (band.end <= age) continue;
    // if we are entering a band higher than current age, then age up
    age = std::max(age, band.start);
    // determine how long we will spend in this band
    const double duration = band.end - age;
    // A zero-rate band consumes time but no threshold.
    if (band.rate > 0) {
      if (threshold <= band.rate * duration)
        // The event will occur within this band. Return the total waiting time.
        // Include time crossed in previous bands, then the remaining time here.
        return age - origin + threshold / band.rate;
      //if the event doesn't happen consume hazard from this band
      threshold -= band.rate * duration;
    }
    age = band.end;
  }
  // if we got here, we are out of bands so the event never happens, return infinity
  return std::numeric_limits<double>::infinity();
}

// Internal deterministic entry point for analytical validation of the sampler.
// [[Rcpp::export]]
NumericVector hazard_wait_cpp(DataFrame rates, NumericVector ages,
                             NumericVector thresholds) {
  if (ages.size() != thresholds.size()) stop("Sampler inputs must have equal length.");
  std::vector<Band> bands = read_bands(rates);
  NumericVector waits(ages.size());
  for (int i = 0; i < ages.size(); ++i) {
    if (!R_finite(ages[i]) || ages[i] < 0 || !R_finite(thresholds[i]) ||
        thresholds[i] <= 0) stop("Sampler ages and thresholds are invalid.");
    waits[i] = waiting_time(bands, ages[i], thresholds[i]);
  }
  return waits;
}

struct Person {
  // Public history fields plus working state. NA died means currently alive;
  // last_birth enforces spacing. male_slot locates this person in available males.
  double id, born, died, mother, father, last_birth;
  bool female;
  int male_slot, group;
};
// person is a zero-based position in people, not the public person_id.
// A version invalidates queued events after a change in union status.
enum class EventType { death, birth, search };
struct Event { int person; EventType type; unsigned long long version; };

double rate_integral(const std::vector<Band>& bands, double lo, double hi) {
  // Integrate the supplied annual rate over [lo, hi), accounting for bands
  // that only partly overlap the interval. Uncovered fertility ages add zero.
  double total = 0;
  for (const Band& band : bands) {
    if (band.end <= lo) continue;
    if (band.start >= hi) break;
    total += band.rate * std::max(0.0, std::min(hi, band.end) - std::max(lo, band.start));
  }
  return total;
}

// Calculate one constant eligible hazard at each original band's midpoint.
// Integrate the unadjusted rates across the preceding birth interval, including
// overlaps with earlier bands and gaps. Event sampling uses the ordinary Band
// representation and shared constant-hazard waiting_time() function.
std::vector<Band> adjusted_fertility(const std::vector<Band>& raw,
                                    double birth_interval) {
  for (const Band& band : raw) {
    if (!std::isfinite(band.end)) stop("Fertility age bands must have finite ends.");
  }
  if (birth_interval == 0) return raw;
  std::vector<Band> adjusted;
  for (const Band& band : raw) {
    checkUserInterrupt();
    const double midpoint = band.start + (band.end - band.start) / 2;
    const double eligible = 1 - rate_integral(raw,
      std::max(0.0, midpoint - birth_interval), midpoint);
    if (!(eligible > 0))
      stop("Birth interval and fertility rates imply nonpositive midpoint eligibility.");
    const double rate = band.rate / eligible;
    if (!std::isfinite(rate)) stop("Adjusted fertility rate cannot be represented.");
    adjusted.push_back({band.start, band.end, rate});
  }
  return adjusted;
}

// Internal deterministic entry point to inspect adjustment in analytical tests.
// [[Rcpp::export]]
DataFrame adjusted_fertility_cpp(DataFrame rates, double birth_interval) {
  if (!R_finite(birth_interval) || birth_interval < 0)
    stop("Birth interval must be finite and nonnegative.");
  const auto bands = adjusted_fertility(read_bands(rates), birth_interval);
  NumericVector starts(bands.size()), ends(bands.size()), values(bands.size());
  for (size_t i = 0; i < bands.size(); ++i) {
    starts[i] = bands[i].start;
    ends[i] = bands[i].end;
    values[i] = bands[i].rate;
  }
  return DataFrame::create(_["age_start"] = starts, _["age_end"] = ends,
                          _["rate"] = values);
}

// Partial Fisher-Yates using reusable positions, restored after each sample.
// This avoids allocating a hash node for each candidate at every search.
// Workspace positions are always identity between calls, including when the
// living-male pool grows or shrinks. They index males, not fixed person IDs.
std::vector<int> sample_male_positions(int count, int requested,
                                     std::vector<int>& workspace) {
  const int take = requested == 0 ? count : std::min(count, requested);
  std::vector<int> sample;
  sample.reserve(take);
  if (take == count) {
    for (int j = 0; j < count; ++j) sample.push_back(j);
    return sample;
  }
  while (workspace.size() < static_cast<size_t>(count))
    workspace.push_back(static_cast<int>(workspace.size()));
  std::vector<int> swaps;
  swaps.reserve(take);
  for (int j = 0; j < take; ++j) {
    if (j % 1024 == 0) checkUserInterrupt();
    const int selected = j + std::min(count - j - 1,
      static_cast<int>(R::runif(0, 1) * (count - j)));
    std::swap(workspace[j], workspace[selected]);
    sample.push_back(workspace[j]);
    swaps.push_back(selected);
  }
  // Undo in reverse so the next search starts with identity positions again.
  for (int j = take - 1; j >= 0; --j)
    std::swap(workspace[j], workspace[swaps[j]]);
  return sample;
}

// Evaluate candidate utilities without rounding age differences. Common linear
// and raw quadratic formulas use an exact native path; other formulas call the
// R model-matrix evaluator with candidate data only, never mutable engine state.
NumericVector matching_utilities(const NumericVector& gaps, const List& matching,
                                 int female_group = -1,
                                 IntegerVector male_groups = IntegerVector(),
                                 IntegerVector female_groups = IntegerVector()) {
  NumericVector quadratic = matching["quadratic"];
  NumericVector utility(gaps.size());
  if (quadratic.size() == 2) {
    const bool has_squared = as<bool>(matching["has_squared"]);
    for (R_xlen_t j = 0; j < gaps.size(); ++j) {
      const double squared = has_squared ? gaps[j] * gaps[j] : 0;
      if (!R_finite(squared)) stop("Matching formula must produce finite values for every candidate.");
      utility[j] = quadratic[0] * gaps[j] + quadratic[1] * squared;
    }
  } else {
    Function evaluate = matching["evaluate"];
    if (as<bool>(matching["group_formula"])) {
      if ((female_group < 0 && female_groups.size() != gaps.size()) || male_groups.size() != gaps.size())
        stop("Matching group predictors require group data.");
      if (female_groups.size())
        utility = as<NumericVector>(evaluate(gaps, female_groups, male_groups));
      else utility = as<NumericVector>(evaluate(gaps, female_group, male_groups));
    } else utility = as<NumericVector>(evaluate(gaps));
  }
  if (utility.size() != gaps.size()) stop("Matching formula returned the wrong number of utilities.");
  NumericMatrix offsets = matching["pair_offsets"];
  if (offsets.nrow()) {
    if ((female_groups.size() != gaps.size() && (female_group < 0 || female_group >= offsets.nrow())) || male_groups.size() != gaps.size())
      stop("Pair effects require group data.");
    for (R_xlen_t j = 0; j < utility.size(); ++j) {
      if (male_groups[j] < 0 || male_groups[j] >= offsets.ncol()) stop("Invalid matching group index.");
      const int fg = female_groups.size() ? female_groups[j] : female_group;
      if (fg < 0 || fg >= offsets.nrow()) stop("Invalid matching group index.");
      utility[j] += offsets(fg, male_groups[j]);
    }
  }
  for (R_xlen_t j = 0; j < utility.size(); ++j) {
    if (!R_finite(gaps[j]) || !R_finite(utility[j]))
      stop("Matching ages and utilities must be finite.");
  }
  return utility;
}

NumericVector matching_probabilities(const NumericVector& utility) {
  if (utility.size() == 0) return NumericVector();
  const double maximum = *std::max_element(utility.begin(), utility.end());
  NumericVector probabilities(utility.size());
  double total = 0;
  for (R_xlen_t j = 0; j < utility.size(); ++j) {
    probabilities[j] = std::exp(utility[j] - maximum);
    total += probabilities[j];
  }
  for (double& probability : probabilities) probability /= total;
  return probabilities;
}

// One uniform draw selects a candidate from the softmax cumulative distribution.
int draw_partner_position(const NumericVector& probabilities) {
  if (probabilities.size() == 0) return -1;
  const double draw = R::runif(0, 1);
  double cumulative = 0;
  int last_positive = -1;
  for (R_xlen_t j = 0; j < probabilities.size(); ++j) {
    if (probabilities[j] > 0) last_positive = static_cast<int>(j);
    cumulative += probabilities[j];
    if (draw < cumulative) return static_cast<int>(j);
  }
  // Roundoff must not select a zero-weight tail candidate.
  return last_positive;
}

// Local hypothetical markets never mutate population or event state.
struct TwoSidedResult { int man = -1, steps = 0, chooser = -1; double scored = 0; };

TwoSidedResult finish_two_sided_market(int nw, int nm,
    const std::vector<double>& column_max, const std::vector<double>& column_sum,
    const std::vector<double>& ego_score, TwoSidedResult result) {
  const double maximum = *std::max_element(ego_score.begin(), ego_score.end());
  if (maximum == R_NegInf) return result;
  double row_sum = 0;
  for (double value : ego_score) row_sum += std::exp(value - maximum);
  std::vector<double> weights(nm);
  double total = 0;
  for (int j = 0; j < nm; ++j) {
    if (ego_score[j] == R_NegInf) continue;
    const double a = std::exp(ego_score[j] - maximum) / row_sum;
    const double b = std::exp(ego_score[j] - column_max[j]) / column_sum[j];
    // Equal pool sizes in production; factors also support rectangular diagnostics.
    weights[j] = a / nw + b / nm;
    total += weights[j];
  }
  const double draw = R::runif(0, 1) * total;
  double cumulative = 0;
  for (int j = 0; j < nm; ++j) {
    if (weights[j] > 0) result.man = j;
    cumulative += weights[j];
    if (draw < cumulative) break;
  }
  result.steps = 1; // One integrated draw; no initiator is sampled.
  return result;
}

template <typename Score>
TwoSidedResult draw_two_sided_market(int nw, int nm, int ego, Score score) {
  TwoSidedResult result;
  if (!nw || !nm) return result;
  std::vector<int> women(nw), men(nm);
  for (int i = 0; i < nw; ++i) women[i] = i;
  for (int j = 0; j < nm; ++j) men[j] = j;
  std::vector<double> column_max(nm, R_NegInf), column_sum(nm, 0), ego_score(nm);
  for (int i = 0; i < nw; ++i) {
    if (i % 64 == 0) checkUserInterrupt();
    const NumericVector values = score(women, men, true, i);
    if (values.size() != nm) stop("Temporary matching returned the wrong number of scores.");
    for (int j = 0; j < nm; ++j) {
      const double value = values[j];
      ++result.scored;
      if (i == ego) ego_score[j] = value;
      if (value == R_NegInf) continue;
      if (!R_finite(value)) stop("Temporary matching utilities must be finite or excluded.");
      if (value > column_max[j]) {
        column_sum[j] = column_sum[j] * std::exp(column_max[j] - value) + 1;
        column_max[j] = value;
      } else column_sum[j] += std::exp(value - column_max[j]);
    }
  }
  return finish_two_sided_market(nw, nm, column_max, column_sum, ego_score, result);
}

// Evaluate each fixed female row once; reuse its scores in both normalizations.
class TemporaryScores {
  NumericVector female_born, male_born, quadratic;
  IntegerVector female_groups, male_groups;
  NumericMatrix offsets;
  LogicalMatrix blocked;
  List matching;
  bool native, squared, prevent_all;
  bool allowed(int i, int j) const {
    const int fg = female_groups[i], mg = male_groups[j];
    return !(prevent_all && fg != mg) && !(blocked.nrow() && blocked(fg, mg));
  }
 public:
  TemporaryScores(NumericVector fb, NumericVector mb, IntegerVector fg,
                  IntegerVector mg, const List& model) :
      female_born(fb), male_born(mb), quadratic(model["quadratic"]),
      female_groups(fg), male_groups(mg), offsets(model["pair_offsets"]),
      blocked(model["blocked_pairs"]), matching(model), native(quadratic.size() == 2),
      squared(as<bool>(model["has_squared"])),
      prevent_all(as<bool>(model["prevent_exogamy_all"])) {
    if (fg.size() != fb.size() || mg.size() != mb.size()) stop("Temporary market group lengths differ.");
    if (as<bool>(model["needs_groups"])) {
      const int groups = std::max(offsets.nrow(), blocked.nrow());
      for (int group : fg) if (group < 0 || (groups && group >= groups)) stop("Invalid female group.");
      for (int group : mg) if (group < 0 || (groups && group >= groups)) stop("Invalid male group.");
    }
  }
  bool is_native() const { return native; }

  // Batched matrix columns: no R allocations or Rcpp element proxies per row.
  // Keep raw gap arithmetic, rather than expanding calendar-time squares.
  TwoSidedResult draw_native(int ego = 0) const {
    const int nw = female_born.size(), nm = male_born.size();
    TwoSidedResult result;
    if (!nw || !nm) return result;
    const double* fb = REAL(female_born);
    const double* mb = REAL(male_born);
    const int* fg = INTEGER(female_groups);
    const int* mg = INTEGER(male_groups);
    const int ng = offsets.nrow(), nb = blocked.nrow();
    const double* pair = ng ? REAL(offsets) : nullptr;
    const int* excluded = nb ? LOGICAL(blocked) : nullptr;
    const double linear = quadratic[0], quadratic_term = quadratic[1];
    std::vector<double> column_max(nm, R_NegInf), column_sum(nm, 0), ego_score(nm), values(nw);
    for (int j = 0; j < nm; ++j) {
      if (j % 64 == 0) checkUserInterrupt();
      const double male_birth = mb[j];
      const int male_group = mg[j];
      double maximum = R_NegInf;
      for (int i = 0; i < nw; ++i) {
        if ((prevent_all && fg[i] != male_group) ||
            (excluded && excluded[fg[i] + nb * male_group])) {
          values[i] = R_NegInf;
          continue;
        }
        const double gap = fb[i] - male_birth;
        const double square = squared ? gap * gap : 0;
        if (!R_finite(gap) || !R_finite(square)) stop("Matching ages and utilities must be finite.");
        double value = linear * gap + quadratic_term * square;
        if (pair) value += pair[fg[i] + ng * male_group];
        if (!R_finite(value)) stop("Matching ages and utilities must be finite.");
        values[i] = value;
        maximum = std::max(maximum, value);
      }
      column_max[j] = maximum;
      ego_score[j] = values[ego];
      if (maximum != R_NegInf) {
        double sum = 0;
        for (double value : values) sum += std::exp(value - maximum);
        column_sum[j] = sum;
      }
    }
    result.scored = static_cast<double>(nw) * nm;
    return finish_two_sided_market(nw, nm, column_max, column_sum, ego_score, result);
  }

  NumericVector choice(const std::vector<int>& women, const std::vector<int>& men,
                       bool female_chooses, int initiator) const {
    const int n = female_chooses ? men.size() : women.size();
    NumericVector output(n, R_NegInf);
    std::vector<int> positions;
    for (int k = 0; k < n; ++k) {
      const int i = female_chooses ? women[initiator] : women[k];
      const int j = female_chooses ? men[k] : men[initiator];
      if (!allowed(i, j)) continue;
      if (!native) { positions.push_back(k); continue; }
      const double gap = female_born[i] - male_born[j];
      const double square = squared ? gap * gap : 0;
      if (!R_finite(gap) || !R_finite(square)) stop("Matching ages and utilities must be finite.");
      double value = quadratic[0] * gap + quadratic[1] * square;
      if (offsets.nrow()) value += offsets(female_groups[i], male_groups[j]);
      if (!R_finite(value)) stop("Matching ages and utilities must be finite.");
      output[k] = value;
    }
    if (!native && !positions.empty()) {
      NumericVector gaps(positions.size());
      IntegerVector fg(positions.size()), mg(positions.size());
      for (size_t k = 0; k < positions.size(); ++k) {
        const int i = female_chooses ? women[initiator] : women[positions[k]];
        const int j = female_chooses ? men[positions[k]] : men[initiator];
        gaps[k] = female_born[i] - male_born[j]; fg[k] = female_groups[i]; mg[k] = male_groups[j];
      }
      const auto values = matching_utilities(gaps, matching, -1, mg, fg);
      for (size_t k = 0; k < positions.size(); ++k) output[positions[k]] = values[k];
    }
    return output;
  }
};

// Internal deterministic-market entry points use the actual simulation kernel.
// [[Rcpp::export]]
List inspect_two_sided_matrix_cpp(NumericMatrix utility, int ego, int draws) {
  if (ego < 1 || ego > utility.nrow() || draws < 0) stop("Invalid ego or number of draws.");
  for (double value : utility)
    if (!R_finite(value) && value != R_NegInf) stop("Utilities must be finite or -Inf for exclusions.");
  IntegerVector choices(draws), steps(draws), chooser(draws);
  NumericVector scored(draws);
  for (int k = 0; k < draws; ++k) {
    if (k % 1024 == 0) checkUserInterrupt();
    auto result = draw_two_sided_market(utility.nrow(), utility.ncol(), ego - 1,
      [&](const std::vector<int>& women, const std::vector<int>& men, bool female, int initiator) {
        NumericVector values(female ? men.size() : women.size());
        for (R_xlen_t k = 0; k < values.size(); ++k)
          values[k] = female ? utility(women[initiator], men[k]) : utility(women[k], men[initiator]);
        return values;
      });
    choices[k] = result.man < 0 ? NA_INTEGER : result.man + 1;
    steps[k] = result.steps; chooser[k] = result.chooser; scored[k] = result.scored;
  }
  return List::create(_["choices"] = choices, _["steps"] = steps,
                      _["chooser"] = chooser, _["scored_pairs"] = scored);
}

// [[Rcpp::export]]
List inspect_two_sided_matching_cpp(NumericVector female_born, NumericVector male_born,
                                    IntegerVector female_groups, IntegerVector male_groups,
                                    List matching, int draws, bool precompute = false) {
  if (!female_born.size() || draws < 0) stop("Need ego and a nonnegative number of draws.");
  TemporaryScores scores(female_born, male_born, female_groups, male_groups, matching);
  IntegerVector choices(draws), steps(draws), chooser(draws);
  NumericVector scored(draws);
  for (int k = 0; k < draws; ++k) {
    if (k % 1024 == 0) checkUserInterrupt();
    // Rebuild the eager matrix for every draw: benchmarking one fresh market
    // must not amortize its cost over repeated selections from a fixed matrix.
    std::vector<double> grid;
    if (precompute) {
      std::vector<int> women(female_born.size()), men(male_born.size());
      for (int i = 0; i < female_born.size(); ++i) women[i] = i;
      for (int j = 0; j < male_born.size(); ++j) men[j] = j;
      grid.reserve(static_cast<size_t>(women.size()) * men.size());
      for (int i = 0; i < female_born.size(); ++i) {
        if (i % 64 == 0) checkUserInterrupt();
        const auto row = scores.choice(women, men, true, i);
        grid.insert(grid.end(), row.begin(), row.end());
      }
    }
    auto result = (!precompute && scores.is_native()) ? scores.draw_native() :
      draw_two_sided_market(female_born.size(), male_born.size(), 0,
      [&](const std::vector<int>& women, const std::vector<int>& men, bool female, int initiator) {
        if (!precompute) return scores.choice(women, men, female, initiator);
        NumericVector values(female ? men.size() : women.size());
        for (R_xlen_t j = 0; j < values.size(); ++j) {
          const int wi = female ? women[initiator] : women[j];
          const int mj = female ? men[j] : men[initiator];
          values[j] = grid[static_cast<size_t>(wi) * male_born.size() + mj];
        }
        return values;
      });
    choices[k] = result.man < 0 ? NA_INTEGER : result.man + 1;
    steps[k] = result.steps; chooser[k] = result.chooser; scored[k] = result.scored;
  }
  return List::create(_["choices"] = choices, _["steps"] = steps,
                      _["chooser"] = chooser, _["scored_pairs"] = scored);
}

// Apply search eligibility after sampling; retain positions in the male pool.
// Further pair restrictions can be added here without changing sample size.
std::vector<int> screen_male_candidates(const std::vector<int>& males,
                                       const std::vector<int>& sample,
                                       const std::vector<Person>& people,
                                       double now, double min_male_age, int woman,
                                       hammel::Kinship& kinship) {
  std::vector<int> eligible;
  eligible.reserve(sample.size());
  const bool check_kinship = kinship.generations() > 0;
  for (int position : sample) {
    if (now - people[males[position]].born >= min_male_age &&
        (!check_kinship || !kinship.known_related(woman, males[position])))
      eligible.push_back(position);
  }
  return eligible;
}

// Select a partner from sampled positions in males, returning a people position.
// Age(male)-age(female) equals birth_time(female)-birth_time(male).
int select_partner(const std::vector<int>& males, const std::vector<int>& sample,
                   const std::vector<Person>& people, int woman, const List& matching,
                   int* utility_rank = nullptr) {
  if (sample.empty()) return -1;
  NumericVector gaps(sample.size());
  for (size_t j = 0; j < sample.size(); ++j)
    gaps[j] = people[woman].born - people[males[sample[j]]].born;
  IntegerVector male_groups;
  if (as<bool>(matching["needs_groups"])) {
    male_groups = IntegerVector(sample.size());
    for (size_t j = 0; j < sample.size(); ++j)
      male_groups[j] = people[males[sample[j]]].group;
  }
  const NumericVector utility = matching_utilities(gaps, matching, people[woman].group, male_groups);
  const NumericVector probabilities = matching_probabilities(utility);
  const int chosen = draw_partner_position(probabilities);
  if (utility_rank) {
    *utility_rank = 1;
    for (double score : utility)
      if (score > utility[chosen]) ++*utility_rank;
  }
  return males[sample[chosen]];
}

// Internal analytical/sampling diagnostics use the same scoring and draw helpers.
// [[Rcpp::export]]
List inspect_matching_cpp(NumericVector age_difference, List matching, int draws) {
  if (draws < 0) stop("Number of draws must be nonnegative.");
  const NumericVector utility = matching_utilities(age_difference, matching);
  const NumericVector probabilities = matching_probabilities(utility);
  IntegerVector choices(draws);
  for (int j = 0; j < draws; ++j) {
    if (j % 1024 == 0) checkUserInterrupt();
    const int selected = draw_partner_position(probabilities);
    choices[j] = selected < 0 ? NA_INTEGER : selected + 1;
  }
  return List::create(_["utilities"] = utility, _["probabilities"] = probabilities,
                      _["choices"] = choices);
}

// Analytical diagnostics for group-aware scoring use the same native kernel.
// Group indices refer to the labels supplied to prepare_matching_model().
// [[Rcpp::export]]
List inspect_group_matching_cpp(NumericVector age_difference, List matching,
                                int female_group, IntegerVector male_groups) {
  const NumericVector utility = matching_utilities(age_difference, matching,
                                                  female_group, male_groups);
  return List::create(_["utilities"] = utility,
                      _["probabilities"] = matching_probabilities(utility));
}

// Internal sampling entry point for checking size, uniqueness, and inclusion.
// [[Rcpp::export]]
IntegerVector sample_male_positions_cpp(int count, int requested) {
  if (count < 0 || requested < 0) stop("Sample sizes must be nonnegative.");
  std::vector<int> workspace;
  return wrap(sample_male_positions(count, requested, workspace));
}

struct UnionRecord {
  double id, start, end;
  int wife, husband;
  bool wife_died = false;
};

// [[Rcpp::export]]
List simulate_cpp(DataFrame persons, DataFrame mortality_rates,
                  DataFrame fertility_rates, double start, double years,
                  double birth_interval, DataFrame unions,
                  int sample_size, List matching, double min_male_age,
                  int kinship_generations, bool rate_initiation,
                  DataFrame union_rates, List father_matching,
                  bool cancel_birth_if_no_match, bool diagnostics,
                  NumericMatrix initial_ancestry, List inheritance, double ffab) {
  if (!R_finite(ffab) || ffab < 0 || ffab > 1)
    stop("ffab must be a finite probability in [0, 1].");
  // R performs input validation and seed handling. All engine state below is
  // local to this call; public times and rates remain in years and person-years.
  NumericVector ids = persons["person_id"], births = persons["birth_time"];
  NumericVector deaths = persons["death_time"], moms = persons["mother_id"];
  NumericVector dads = persons["father_id"];
  CharacterVector sexes = persons["sex"];
  // Row-major ancestry shares belong to this simulation. The last component
  // is unknown; zero components means ancestry tracking is absent.
  const size_t ancestry_components = initial_ancestry.ncol();
  if (initial_ancestry.nrow() != ids.size()) stop("Ancestry rows must match persons.");
  std::vector<double> ancestry(static_cast<size_t>(ids.size()) * ancestry_components);
  for (int i = 0; i < ids.size(); ++i) {
    for (size_t g = 0; g < ancestry_components; ++g)
      ancestry[static_cast<size_t>(i) * ancestry_components + g] = initial_ancestry(i, g);
  }
  IntegerVector person_groups(ids.size(), 0);
  if (persons.containsElementNamed(".hammel_group")) person_groups = persons[".hammel_group"];
  // people retains every history, including deaths. males contains positions
  // of living, unpartnered males for candidate sampling. index maps
  // public IDs to positions while reconstructing initial maternal histories.
  std::vector<Person> people;
  std::vector<int> males;
  std::vector<int> sampling_workspace;
  std::map<double, int> index;
  double next_id = 0;
  for (int i = 0; i < ids.size(); ++i) {
    bool female = as<std::string>(sexes[i]) == "f";
    bool alive = NumericVector::is_na(deaths[i]);
    int slot = (!female && alive) ? static_cast<int>(males.size()) : -1;
    if (slot >= 0) males.push_back(i);
    people.push_back({ids[i], births[i], deaths[i], moms[i], dads[i],
                      -std::numeric_limits<double>::infinity(), female, slot, person_groups[i]});
    index[ids[i]] = i;
    next_id = std::max(next_id, ids[i]);
  }
  std::vector<std::array<int, 2>> parent_links(people.size());
  for (int person = 0; person < static_cast<int>(people.size()); ++person) {
    parent_links[person] = {
      NumericVector::is_na(people[person].mother) ? -1 : index.at(people[person].mother),
      NumericVector::is_na(people[person].father) ? -1 : index.at(people[person].father)
    };
  }
  hammel::Kinship kinship(std::move(parent_links), kinship_generations);
  std::vector<UnionRecord> union_history;
  std::vector<int> diagnostic_counts, diagnostic_ranks;
  std::vector<std::vector<int>> person_unions(people.size());
  std::vector<int> current_partner(people.size(), -1);
  // Swap removal keeps pool updates constant-time and slots valid after moves.
  auto remove_male = [&](int man) {
    const int slot = people[man].male_slot;
    if (slot < 0) return; // Partnered men are already absent when they die.
    const int moved = males.back();
    males[slot] = moved;
    people[moved].male_slot = slot;
    males.pop_back();
    people[man].male_slot = -1;
  };
  auto add_male = [&](int man) {
    if (people[man].female || !NumericVector::is_na(people[man].died) ||
        current_partner[man] >= 0 || people[man].male_slot >= 0) return;
    people[man].male_slot = static_cast<int>(males.size());
    males.push_back(man);
  };
  double next_union_id = 0;
  NumericVector uids = unions["union_id"], wives = unions["wife_id"];
  NumericVector husbands = unions["husband_id"], starts = unions["start_time"];
  NumericVector ends = unions["end_time"];
  for (int u = 0; u < uids.size(); ++u) {
    const int woman = index.at(wives[u]), man = index.at(husbands[u]);
    union_history.push_back({uids[u], starts[u], ends[u], woman, man});
    next_union_id = std::max(next_union_id, uids[u]);
    if (!NumericVector::is_na(ends[u])) continue;
    if (!NumericVector::is_na(people[woman].died) ||
        !NumericVector::is_na(people[man].died) ||
        !people[woman].female || people[man].female)
      stop("Unions require living female/male partners in active input unions.");
    if (current_partner[woman] >= 0)
      stop("Union input cannot have multiple ongoing unions for a woman.");
    if (current_partner[man] >= 0)
      stop("Union input cannot have multiple ongoing unions for a man.");
    current_partner[woman] = man;
    current_partner[man] = woman;
    remove_male(man);
    person_unions[woman].push_back(u);
    person_unions[man].push_back(u);
  }
  // On restart, recover the latest recorded maternal birth, including children
  // who have died. A mother with no recorded children has no spacing restriction.
  for (const Person& child : people) {
    if (!NumericVector::is_na(child.mother)) {
      auto found = index.find(child.mother);
      if (found != index.end())
        people[found->second].last_birth = std::max(people[found->second].last_birth, child.born);
    }
  }
  const bool two_sided = as<bool>(matching["two_sided"]);
  std::vector<int> females, female_slots(two_sided ? people.size() : 0, -1);
  std::vector<int> female_sampling_workspace, reduction_workspace;
  std::multimap<double, int> female_age_entries;
  auto remove_female = [&](int woman) {
    if (!two_sided || female_slots[woman] < 0) return;
    const int slot = female_slots[woman], moved = females.back();
    females[slot] = moved; female_slots[moved] = slot;
    females.pop_back(); female_slots[woman] = -1;
  };
  auto add_female = [&](int woman, double now) {
    if (!two_sided || !people[woman].female || !NumericVector::is_na(people[woman].died) ||
        current_partner[woman] >= 0 || female_slots[woman] >= 0) return;
    if (now - people[woman].born < min_male_age) {
      double entry = people[woman].born + min_male_age;
      // Avoid rescheduling at an already-reached rounded birthday when the
      // subtraction-based age check still places her just below the cutoff.
      if (entry <= now) entry = std::nextafter(now, std::numeric_limits<double>::infinity());
      female_age_entries.emplace(entry, woman);
      return;
    }
    female_slots[woman] = females.size(); females.push_back(woman);
  };
  if (two_sided) {
    for (int woman = 0; woman < static_cast<int>(people.size()); ++woman) add_female(woman, start);
  }
  // All father and partner searches share sampling, eligibility and matching.
  // Only form_union changes partner state or removes the selected man.
  auto find_candidate = [&](int woman, double now, const List& choice_model,
                            int* count = nullptr, int* rank = nullptr) {
    if (males.empty()) return -1;
    const auto sample = sample_male_positions(static_cast<int>(males.size()), sample_size,
                                             sampling_workspace);
    auto eligible = screen_male_candidates(males, sample, people, now, min_male_age,
                                           woman, kinship);
    const bool prevent_all = choice_model["prevent_exogamy_all"];
    LogicalMatrix blocked = choice_model["blocked_pairs"];
    if (prevent_all || blocked.nrow()) {
      const int female_group = people[woman].group;
      eligible.erase(std::remove_if(eligible.begin(), eligible.end(), [&](int position) {
        const int male_group = people[males[position]].group;
        return (prevent_all && female_group != male_group) ||
          (blocked.nrow() && blocked(female_group, male_group));
      }), eligible.end());
    }
    if (count) *count = static_cast<int>(eligible.size());
    if (!two_sided || eligible.empty())
      return select_partner(males, eligible, people, woman, choice_model, rank);
    // Promote age-eligible women before sampling, without adding execution steps.
    while (!female_age_entries.empty() && female_age_entries.begin()->first <= now) {
      const auto entry = *female_age_entries.begin();
      female_age_entries.erase(female_age_entries.begin());
      add_female(entry.second, now);
    }
    const int ego_slot = female_slots[woman];
    const int available = static_cast<int>(females.size()) - (ego_slot >= 0 ? 1 : 0);
    const int market_size = std::min(static_cast<int>(eligible.size()), available + 1);
    const auto original_eligible = rank ? eligible : std::vector<int>();
    if (market_size < static_cast<int>(eligible.size())) {
      const auto kept = sample_male_positions(eligible.size(), market_size, reduction_workspace);
      std::vector<int> reduced;
      reduced.reserve(market_size);
      for (int k : kept) reduced.push_back(eligible[k]);
      eligible.swap(reduced);
    }
    std::vector<int> women{woman};
    if (market_size > 1) {
      const auto positions = sample_male_positions(available, market_size - 1, female_sampling_workspace);
      for (int position : positions)
        women.push_back(females[position + (ego_slot >= 0 && position >= ego_slot ? 1 : 0)]);
    }
    NumericVector fb(market_size), mb(market_size);
    IntegerVector fg(market_size), mg(market_size);
    for (int k = 0; k < market_size; ++k) {
      fb[k] = people[women[k]].born; fg[k] = people[women[k]].group;
      mb[k] = people[males[eligible[k]]].born; mg[k] = people[males[eligible[k]]].group;
    }
    TemporaryScores scores(fb, mb, fg, mg, choice_model);
    const auto outcome = scores.is_native() ? scores.draw_native() :
      draw_two_sided_market(market_size, market_size, 0,
      [&](const std::vector<int>& women, const std::vector<int>& men, bool female, int initiator) {
        return scores.choice(women, men, female, initiator);
      });
    if (outcome.man < 0) stop("Ego lost all eligible temporary partners.");
    const int selected = males[eligible[outcome.man]];
    if (rank) {
      NumericVector gaps(original_eligible.size());
      IntegerVector groups(original_eligible.size());
      int selected_position = -1;
      for (size_t j = 0; j < original_eligible.size(); ++j) {
        const int man = males[original_eligible[j]];
        gaps[j] = people[woman].born - people[man].born; groups[j] = people[man].group;
        if (man == selected) selected_position = j;
      }
      const auto utility = matching_utilities(gaps, choice_model, people[woman].group, groups);
      *rank = 1;
      for (double value : utility) if (value > utility[selected_position]) ++*rank;
    }
    return selected;
  };
  auto form_union = [&](int woman, double now) {
    int partner = current_partner[woman];
    if (partner >= 0) return partner;
    int count = NA_INTEGER, rank = NA_INTEGER;
    partner = find_candidate(woman, now, matching,
                             diagnostics ? &count : nullptr,
                             diagnostics ? &rank : nullptr);
    if (partner >= 0) {
      if (next_union_id >= 9007199254740991.0 ||
          union_history.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Union history exceeds supported ID or row limits.");
      const int u = static_cast<int>(union_history.size());
      union_history.push_back({++next_union_id, now, NA_REAL, woman, partner});
      if (diagnostics) {
        diagnostic_counts.push_back(count);
        diagnostic_ranks.push_back(rank);
      }
      person_unions[woman].push_back(u);
      person_unions[partner].push_back(u);
      current_partner[woman] = partner;
      current_partner[partner] = woman;
      remove_male(partner);
      remove_female(woman);
    }
    return partner;
  };
  struct GroupSchedules {
    std::vector<Band> male, female, partnered, unpartnered, search;
  };
  int group_count = std::max(1, as<int>(inheritance["group_count"]));
  const int inheritance_method = inheritance["method"];
  const double passing = inheritance["passing"];
  IntegerVector inheritance_ranks = inheritance["ranks"];
  IntegerVector ancestry_groups = inheritance["ancestry_groups"];
  Function inherit_custom = inheritance["custom"];
  for (int value : person_groups) group_count = std::max(group_count, value + 1);
  std::vector<GroupSchedules> schedules(group_count);
  const bool status_fertility = fertility_rates.containsElementNamed("union_status");
  for (int group = 0; group < group_count; ++group) {
    auto& rates = schedules[group];
    rates.male = read_bands(mortality_rates, "m", "", group);
    rates.female = read_bands(mortality_rates, "f", "", group);
    rates.partnered = adjusted_fertility(read_bands(fertility_rates, "",
      status_fertility ? "partnered" : "", group), birth_interval);
    rates.unpartnered = status_fertility ? adjusted_fertility(
      read_bands(fertility_rates, "", "unpartnered", group), birth_interval) : rates.partnered;
    if (rate_initiation) rates.search = read_bands(union_rates, "", "", group);
  }
  std::vector<unsigned long long> event_versions(people.size(), 0);
  // Continuous dates order events. Equivalent keys retain insertion order.
  std::multimap<double, Event> queue;
  const double end = start + years;
  auto candidate = [&](double waiting, double now) {
    double due = now + waiting;
    if (std::isfinite(waiting) && due <= now)
      stop("Event time cannot advance at this numerical scale.");
    return due;
  };
  auto threshold = [&]() {
    // Use R's RNG, shared with the seed policy in simulate(). Exp(1) is the
    // cumulative-hazard threshold; avoid zero to keep waiting times positive.
    double value;
    do { value = R::rexp(1.0); } while (value <= 0);
    return value;
  };
  auto schedule = [&](int i, double now) {
    const auto version = ++event_versions[i];
    const Person& p = people[i];
    if (!NumericVector::is_na(p.died)) return;
    const auto& rates = schedules[p.group];
    const auto& death_rates = p.female ? rates.female : rates.male;
    double death_time = candidate(waiting_time(death_rates, now - p.born, threshold()), now);
    double birth_time = std::numeric_limits<double>::infinity();
    if (p.female) {
      const auto& fertility = current_partner[i] >= 0 ? rates.partnered : rates.unpartnered;
      double first = std::max(now, p.last_birth + birth_interval);
      if (first < end)
        birth_time = candidate(waiting_time(fertility,
          first - p.born, threshold()), first);
    }
    double search_time = std::numeric_limits<double>::infinity();
    if (rate_initiation && p.female && current_partner[i] < 0)
      search_time = candidate(waiting_time(rates.search, now - p.born, threshold()), now);
    const double due = std::min(death_time, std::min(birth_time, search_time));
    // Exact same-person ties select death, then birth, then search.
    const auto type = death_time <= due ? EventType::death :
      (birth_time <= search_time ? EventType::birth : EventType::search);
    if (due < end) queue.emplace(due, Event{i, type, version});
  };
  // Initialize one next event for each living person. Known age changes are
  // already incorporated, so birthdays alone do not trigger resampling.
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    schedule(i, start);
  }
  double birth_count = 0, death_count = 0;
  unsigned long long processed = 0;
  // Process events in chronological order until no scheduled events remain.
  while (!queue.empty()) {
    if (++processed % 1024 == 0) checkUserInterrupt();
    auto next = queue.begin();
    const double now = next->first;
    const Event event = next->second;
    queue.erase(next);
    const int i = event.person;
    if (event.version != event_versions[i]) continue;
    if (event.type == EventType::death) {
      // Retain the person's history, record death, and stop their scheduling.
      people[i].died = now;
      for (int u : person_unions[i]) {
        UnionRecord& relationship = union_history[u];
        if (!NumericVector::is_na(relationship.end)) continue;
        relationship.end = now;
        relationship.wife_died = i == relationship.wife;
        current_partner[relationship.wife] = -1;
        current_partner[relationship.husband] = -1;
        add_male(relationship.husband);
        add_female(relationship.wife, now);
        // Her old birth candidate is invalid once her husband has died.
        if (rate_initiation && i == relationship.husband)
          schedule(relationship.wife, now);
      }
      if (!people[i].female) remove_male(i);
      else remove_female(i);
      ++death_count;
    } else if (event.type == EventType::search) {
      form_union(i, now);
      schedule(i, now);
    } else {
      double father = NA_REAL;
      int partner = current_partner[i];
      if (!rate_initiation && partner < 0) partner = form_union(i, now);
      // A nonmarital father stays available: parentage alone is not a union.
      const int father_position = partner >= 0 ? partner :
        (rate_initiation ? find_candidate(i, now, father_matching) : -1);
      if (father_position < 0 && cancel_birth_if_no_match) {
        // No child was born: retain last_birth and redraw risks from this instant.
        schedule(i, now);
        continue;
      }
      if (father_position >= 0) father = people[father_position].id;
      // If no sampled men are eligible, preserve the birth with unknown father and no union.
      if (next_id >= 9007199254740991.0 ||
          people.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Population exceeds supported ID or row limits.");
      // IDs increase beyond existing histories. Store newborn parentage and
      // sex; union formation above happens before the newborn exists.
      double mother = people[i].id;
      people[i].last_birth = now;
      // Default ffab is 0.4886 - its the signature of Berkeley training!
      bool female = R::runif(0, 1) < ffab;
      int child = static_cast<int>(people.size());
      int slot = female ? -1 : static_cast<int>(males.size());
      people.push_back({++next_id, now, NA_REAL, mother, father,
                        -std::numeric_limits<double>::infinity(), female, slot, people[i].group});
      // Append ancestry before scheduling the child. Missing fathers contribute
      // their entire half to unknown; supplied parental unknown shares transmit
      // through the same deterministic averaging rule.
      const size_t child_offset = static_cast<size_t>(child) * ancestry_components;
      ancestry.resize(child_offset + ancestry_components);
      for (size_t g = 0; g < ancestry_components; ++g) {
        const double paternal = father_position >= 0 ?
          ancestry[static_cast<size_t>(father_position) * ancestry_components + g] :
          (g + 1 == ancestry_components ? 1.0 : 0.0);
        ancestry[child_offset + g] = 0.5 *
          ancestry[static_cast<size_t>(i) * ancestry_components + g] + 0.5 * paternal;
      }
      // Sex and ancestry are now available; assign identity before any newborn
      // hazard or matching lookup. Existing people's identities never change.
      const int maternal_group = people[i].group;
      const int paternal_group = father_position < 0 ? -1 : people[father_position].group;
      int child_group = maternal_group;
      if (inheritance_method == 1 && paternal_group >= 0) child_group = paternal_group;
      else if (inheritance_method == 2 && paternal_group >= 0 && paternal_group != maternal_group)
        child_group = R::runif(0, 1) < 0.5 ? maternal_group : paternal_group;
      else if (inheritance_method == 3 && !female && paternal_group >= 0) child_group = paternal_group;
      else if (inheritance_method == 4 && female && paternal_group >= 0) child_group = paternal_group;
      else if (inheritance_method == 5 && ancestry_components > 0 &&
               ancestry[child_offset + ancestry_components - 1] < 1.0) {
        double total = 0;
        int positive = 0, sole_group = maternal_group;
        for (size_t g = 0; g + 1 < ancestry_components; ++g) {
          const double share = ancestry[child_offset + g];
          if (share > 0) {
            if (ancestry_groups[g] < 0) stop("Positive ancestry has no declared identity group.");
            total += share; ++positive; sole_group = ancestry_groups[g];
          }
        }
        if (positive == 1) child_group = sole_group;
        else if (positive > 1) {
          const double draw = R::runif(0, total);
          double cumulative = 0;
          child_group = sole_group; // Protect the final cell against rounding.
          for (size_t g = 0; g + 1 < ancestry_components; ++g) {
            cumulative += ancestry[child_offset + g];
            if (draw < cumulative) { child_group = ancestry_groups[g]; break; }
          }
        }
      } else if (inheritance_method == 6 && paternal_group >= 0 && paternal_group != maternal_group) {
        const int maternal_rank = inheritance_ranks[maternal_group];
        const int paternal_rank = inheritance_ranks[paternal_group];
        if (maternal_rank >= 0 && paternal_rank >= 0) {
          const bool reverse = passing == 1 || (passing > 0 && R::runif(0, 1) < passing);
          const bool maternal_first = maternal_rank < paternal_rank;
          child_group = (maternal_first != reverse) ? maternal_group : paternal_group;
        }
      } else if (inheritance_method == 7) {
        NumericVector maternal(ancestry_components), paternal(ancestry_components), child_shares(ancestry_components);
        for (size_t g = 0; g < ancestry_components; ++g) {
          maternal[g] = ancestry[static_cast<size_t>(i) * ancestry_components + g];
          if (father_position >= 0)
            paternal[g] = ancestry[static_cast<size_t>(father_position) * ancestry_components + g];
          child_shares[g] = ancestry[child_offset + g];
        }
        // R's random functions must see all preceding native draws, and native
        // scheduling must resume from the state left by the R callback.
        PutRNGstate();
        try {
          child_group = as<int>(inherit_custom(maternal_group, paternal_group,
            female ? "f" : "m", maternal, paternal, child_shares, partner >= 0));
        } catch (...) {
          GetRNGstate();
          throw;
        }
        GetRNGstate();
      }
      if (child_group < 0 || child_group >= group_count) stop("Invalid newborn group index.");
      people[child].group = child_group;
      kinship.append(i, father_position);
      if (!female) males.push_back(child);
      person_unions.emplace_back();
      current_partner.push_back(-1);
      if (two_sided) { female_slots.push_back(-1); add_female(child, now); }
      event_versions.push_back(0);
      ++birth_count;
      // Both newborn and mother need fresh candidates. The mother now has a
      // new last_birth; spacing delays her next birth, but not a possible death.
      schedule(child, now);
      schedule(i, now);
    }
  }
  // Materialize ordinary R columns from the C++ records. The R wrapper restores
  // custom columns, closes existing unions affected by deaths, and assembles the
  // final population and run metadata. Pending events are not saved as a checkpoint.
  NumericVector out_id(people.size()), out_birth(people.size()), out_death(people.size());
  NumericVector out_mother(people.size()), out_father(people.size());
  CharacterVector out_sex(people.size());
  IntegerVector out_group(people.size());
  for (size_t i = 0; i < people.size(); ++i) {
    const Person& p = people[i];
    out_id[i] = p.id; out_birth[i] = p.born; out_death[i] = p.died;
    out_mother[i] = p.mother; out_father[i] = p.father;
    out_sex[i] = p.female ? "f" : "m";
    out_group[i] = p.group;
  }
  NumericMatrix out_ancestry(people.size(), ancestry_components);
  for (size_t i = 0; i < people.size(); ++i) {
    for (size_t g = 0; g < ancestry_components; ++g)
      out_ancestry(i, g) = ancestry[i * ancestry_components + g];
  }
  NumericVector out_uid(union_history.size()), out_wife(union_history.size());
  NumericVector out_husband(union_history.size()), out_start(union_history.size());
  NumericVector out_end(union_history.size());
  CharacterVector out_reason(union_history.size());
  CharacterVector input_reasons = unions["end_reason"];
  for (size_t u = 0; u < union_history.size(); ++u) {
    const UnionRecord& relationship = union_history[u];
    out_uid[u] = relationship.id;
    out_wife[u] = people[relationship.wife].id;
    out_husband[u] = people[relationship.husband].id;
    out_start[u] = relationship.start;
    out_end[u] = relationship.end;
    if (NumericVector::is_na(relationship.end)) out_reason[u] = NA_STRING;
    else if (u < static_cast<size_t>(unions.nrows())) {
      NumericVector input_ends = unions["end_time"];
      if (!NumericVector::is_na(input_ends[u])) out_reason[u] = input_reasons[u];
      else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
    } else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
  }
  DataFrame out_unions = DataFrame::create(_["union_id"] = out_uid,
    _["wife_id"] = out_wife, _["husband_id"] = out_husband,
    _["start_time"] = out_start, _["end_time"] = out_end, _["end_reason"] = out_reason);
  if (diagnostics) {
    IntegerVector counts(union_history.size(), NA_INTEGER), ranks(union_history.size(), NA_INTEGER);
    for (int u = 0; u < unions.nrows(); ++u) {
      if (unions.containsElementNamed(".n_eligible"))
        counts[u] = as<IntegerVector>(unions[".n_eligible"])[u];
      if (unions.containsElementNamed(".utility_rank"))
        ranks[u] = as<IntegerVector>(unions[".utility_rank"])[u];
    }
    for (size_t j = 0; j < diagnostic_counts.size(); ++j) {
      counts[unions.nrows() + j] = diagnostic_counts[j];
      ranks[unions.nrows() + j] = diagnostic_ranks[j];
    }
    out_unions[".n_eligible"] = counts;
    out_unions[".utility_rank"] = ranks;
    out_unions.attr("class") = "data.frame";
    out_unions.attr("row.names") = IntegerVector::create(NA_INTEGER, -static_cast<int>(union_history.size()));
  }
  return List::create(_["persons"] = DataFrame::create(
    _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
    _["death_time"] = out_death, _["mother_id"] = out_mother, _["father_id"] = out_father),
    _["groups"] = out_group, _["unions"] = out_unions, _["ancestry"] = out_ancestry,
    _["births"] = birth_count, _["deaths"] = death_count);
}
