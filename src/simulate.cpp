#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <functional>
#include <memory>
#include <set>
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
                             const std::string& status = "", int group = -1,
                             const std::string& previous = "", bool duration = false, double parity = -1) {
  NumericVector starts = rates[duration ? "dur_start" : "age_start"],
                ends = rates[duration ? "dur_end" : "age_end"];
  NumericVector values = rates["rate"];
  CharacterVector sexes;
  if (!sex.empty()) sexes = rates["sex"];
  CharacterVector statuses;
  if (!status.empty()) statuses = rates["union_type"];
  CharacterVector previous_types;
  if (!previous.empty()) previous_types = rates["previous_union"];
  IntegerVector groups;
  const bool grouped = rates.containsElementNamed(".hammel_group");
  if (grouped) groups = rates[".hammel_group"];
  NumericVector parity_start, parity_end;
  const bool select_parity = parity >= 0 && rates.containsElementNamed("parity_start");
  if (select_parity) { parity_start = rates["parity_start"]; parity_end = rates["parity_end"]; }
  std::vector<Band> bands;
  for (int i = 0; i < starts.size(); ++i) {
    if (select_parity && !(parity_start[i] <= parity && parity < parity_end[i])) continue;
    if ((sex.empty() || as<std::string>(sexes[i]) == sex) &&
        (status.empty() || as<std::string>(statuses[i]) == status) &&
        (!grouped || groups[i] == group) &&
        (previous.empty() || as<std::string>(previous_types[i]) == previous))
      bands.push_back({starts[i], ends[i], values[i]});
  }
  std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) {
    return a.start < b.start;
  });
  return bands;
}

// Select active flat-table rows at a calendar boundary. The public table is
// never changed. R validation checks every active slice before execution.
DataFrame calendar_slice(DataFrame rates, double now) {
  if (!rates.containsElementNamed("time_start")) return rates;
  NumericVector lo = rates["time_start"], hi = rates["time_end"];
  std::vector<int> rows;
  for (int i = 0; i < rates.nrows(); ++i)
    if (lo[i] <= now && now < hi[i]) rows.push_back(i + 1);
  Function subset = Environment::base_env()["["];
  return as<DataFrame>(subset(rates, wrap(rows), R_MissingArg, _["drop"] = false));
}

bool calendar_changes_at(DataFrame rates, double now) {
  if (!rates.containsElementNamed("time_start")) return false;
  NumericVector lo = rates["time_start"], hi = rates["time_end"];
  for (int i = 0; i < rates.nrows(); ++i)
    if (lo[i] == now || hi[i] == now) return true;
  return false;
}

// R validation guarantees one clock per previous category across all groups.
bool search_uses_duration(DataFrame rates, const std::string& previous) {
  if (!rates.containsElementNamed("dur_start") || previous.empty()) return false;
  CharacterVector labels = rates["previous_union"];
  NumericVector ages = rates["age_start"];
  for (int i = 0; i < labels.size(); ++i)
    if (as<std::string>(labels[i]) == previous) return NumericVector::is_na(ages[i]);
  return false;
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
  // last_birth enforces spacing.
  double id, born, died, mother, father, last_birth;
  bool female;
  int group;
};
// person is a zero-based position in people, not the public person_id.
// A version invalidates queued events after a change in union status.
enum class EventType { death, birth, search, dissolution, conversion, transition };
struct Event { int person; EventType type; unsigned long long version; };

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

// Internal analytical diagnostics use the same native/R scoring kernel.
// [[Rcpp::export]]
NumericVector inspect_matching_cpp(NumericVector age_difference, List matching) {
  return matching_utilities(age_difference, matching);
}

// Group indices refer to the labels supplied to prepare_matching_model().
// [[Rcpp::export]]
NumericVector inspect_group_matching_cpp(NumericVector age_difference, List matching,
                                         int female_group, IntegerVector male_groups) {
  return matching_utilities(age_difference, matching, female_group, male_groups);
}

// ---- weighted_choice: calibrated male-group weights (D066) -----------------
// Counts of unpartnered people by group and birth-time bin change only at
// events; ages are evaluated at the search time. Bin width is a numerical
// tolerance, not a demographic parameter. Candidate selection uses exact ages.
constexpr double kWeightBinWidth = 0.25;
constexpr double kZeroTargetLogWeight = -1000; // effectively zero weight
constexpr double kMaxLogWeight = 50;           // cap for unattainable targets
// Incremental pool sums: entries below this floor are recomputed exactly at a
// search, and all sums are rebuilt after this many updates, bounding drift.
constexpr double kPoolSumFloor = 1e-6;
constexpr unsigned long long kPoolSumRefresh = 100000;

double band_rate(const std::vector<Band>& bands, double age) {
  for (const Band& band : bands)
    if (age >= band.start && age < band.end) return band.rate;
  return 0;
}

// Pair utilities on the fixed grid of bin-centre lags. Lag (female bin - male
// bin) * width is male age minus female age, so the grid never changes during
// a run and custom R formulas are evaluated once. Excluded pairs are -Inf.
struct LagTable {
  int groups = 0, span = 0;
  std::vector<double> utility, weight; // weight = exp(utility - global maximum)
  size_t index(int female, int male, int lag) const {
    return (static_cast<size_t>(female) * groups + male) * span + lag;
  }
};

LagTable build_lag_table(const List& model, int groups, int bins, double width) {
  LagTable table;
  table.groups = groups;
  table.span = 2 * bins - 1;
  const size_t n = static_cast<size_t>(groups) * groups * table.span;
  NumericVector gaps(n);
  IntegerVector female(n), male(n);
  for (int f = 0; f < groups; ++f)
    for (int m = 0; m < groups; ++m)
      for (int lag = 0; lag < table.span; ++lag) {
        const size_t i = table.index(f, m, lag);
        gaps[i] = (lag - (bins - 1)) * width;
        female[i] = f;
        male[i] = m;
      }
  const NumericVector utility = matching_utilities(gaps, model, -1, male, female);
  const bool prevent_all = as<bool>(model["prevent_exogamy_all"]);
  LogicalMatrix blocked = model["blocked_pairs"];
  table.utility.assign(utility.begin(), utility.end());
  double maximum = R_NegInf;
  for (int f = 0; f < groups; ++f)
    for (int m = 0; m < groups; ++m) {
      const bool excluded = (prevent_all && f != m) || (blocked.nrow() && blocked(f, m));
      for (int lag = 0; lag < table.span; ++lag) {
        double& value = table.utility[table.index(f, m, lag)];
        if (excluded) value = R_NegInf;
        else maximum = std::max(maximum, value);
      }
    }
  table.weight.resize(n);
  for (size_t i = 0; i < n; ++i)
    table.weight[i] = table.utility[i] == R_NegInf ? 0 : std::exp(table.utility[i] - maximum);
  return table;
}

// Solve a small dense system with partial pivoting; false if singular.
bool solve_linear(std::vector<double> a, std::vector<double> b, int n,
                  std::vector<double>& x) {
  for (int c = 0; c < n; ++c) {
    int pivot = c;
    for (int r = c + 1; r < n; ++r)
      if (std::abs(a[r * n + c]) > std::abs(a[pivot * n + c])) pivot = r;
    if (!(std::abs(a[pivot * n + c]) > 1e-300)) return false;
    if (pivot != c) {
      for (int k = 0; k < n; ++k) std::swap(a[c * n + k], a[pivot * n + k]);
      std::swap(b[c], b[pivot]);
    }
    for (int r = c + 1; r < n; ++r) {
      const double factor = a[r * n + c] / a[c * n + c];
      for (int k = c; k < n; ++k) a[r * n + k] -= factor * a[c * n + k];
      b[r] -= factor * b[c];
    }
  }
  x.assign(n, 0);
  for (int r = n - 1; r >= 0; --r) {
    double sum = b[r];
    for (int k = r + 1; k < n; ++k) sum -= a[r * n + k] * x[k];
    x[r] = sum / a[r * n + r];
  }
  return true;
}

struct WeightSolution {
  std::vector<double> log_weights, target, flow;
  bool infeasible = false;
};

class WeightCalibrator {
  int groups_, bins_;
  long long first_;
  double width_, cutoff_;
  std::vector<double> men_, women_; // [group * bins + bin]
  std::vector<std::array<double, 3>> intensities_; // marriage, fertility
  std::vector<std::array<int, 3>> positive_contributors_;
  // Per registered utility table: sums[(female group * G + male group) * bins
  // + female bin] of count * weight over included male bins. Male bins
  // [0, included_) have reached the age cutoff; older bins have lower indices.
  std::vector<const LagTable*> tables_;
  std::vector<std::vector<double>> sums_;
  int included_ = 0;
  unsigned long long updates_ = 0;
  void spread(size_t t, int male_group, int male_bin, double delta) {
    const LagTable& table = *tables_[t];
    std::vector<double>& sums = sums_[t];
    for (int f = 0; f < groups_; ++f) {
      // Lags for consecutive female bins are contiguous in the table.
      const double* weight = &table.weight[table.index(f, male_group, bins_ - 1 - male_bin)];
      double* out = &sums[(static_cast<size_t>(f) * groups_ + male_group) * bins_];
      for (int k = 0; k < bins_; ++k) out[k] += delta * weight[k];
    }
  }
  void advance(double now) {
    while (included_ < bins_ && now - centre(included_) >= cutoff_) {
      for (int g = 0; g < groups_; ++g) {
        const double count = men_[static_cast<size_t>(g) * bins_ + included_];
        if (count != 0) for (size_t t = 0; t < tables_.size(); ++t) spread(t, g, included_, count);
      }
      ++included_;
    }
  }
  double centre(int k) const { return (static_cast<double>(first_ + k) + 0.5) * width_; }
  size_t cell(int group, double born) const {
    const long long k = static_cast<long long>(std::floor(born / width_)) - first_;
    if (k < 0 || k >= bins_ || group < 0 || group >= groups_)
      stop("Person is outside the weighted-choice calibration range.");
    return static_cast<size_t>(group) * bins_ + static_cast<size_t>(k);
  }

 public:
  WeightCalibrator(double min_birth, double max_birth, int groups, double width, double cutoff)
      : groups_(groups), width_(width), cutoff_(cutoff) {
    if (!(width > 0) || !R_finite(min_birth) || !R_finite(max_birth) || max_birth < min_birth)
      stop("Invalid weighted-choice calibration range.");
    first_ = static_cast<long long>(std::floor(min_birth / width));
    const long long last = static_cast<long long>(std::floor(max_birth / width));
    if (last - first_ + 1 > 50000000LL / std::max(1, groups))
      stop("Weighted-choice calibration range is too large.");
    bins_ = static_cast<int>(last - first_ + 1);
    men_.assign(static_cast<size_t>(groups) * bins_, 0);
    women_ = men_;
    intensities_.resize(men_.size(), {0, 0});
    positive_contributors_.resize(men_.size(), {0, 0});
  }
  int bins() const { return bins_; }
  // Register a utility table before counts are added; returns its index.
  int add_table(const LagTable* table) {
    tables_.push_back(table);
    sums_.emplace_back(static_cast<size_t>(groups_) * groups_ * bins_, 0.0);
    return static_cast<int>(tables_.size()) - 1;
  }
  // Replace a role's preference table; refresh() rebuilds sums from the
  // current pool after all roles have switched. Table storage must stay alive.
  void set_table(int id, const LagTable* table) { tables_.at(id) = table; }
  // Rebuild every pool sum from the counts, discarding accumulated roundoff.
  void refresh() {
    for (size_t t = 0; t < tables_.size(); ++t) {
      std::fill(sums_[t].begin(), sums_[t].end(), 0.0);
      for (int k = 0; k < included_; ++k)
        for (int g = 0; g < groups_; ++g) {
          const double count = men_[static_cast<size_t>(g) * bins_ + k];
          if (count != 0) spread(t, g, k, count);
        }
    }
    updates_ = 0;
  }
  double max_sum_difference(const WeightCalibrator& other, size_t t) const {
    double difference = 0;
    for (size_t i = 0; i < sums_[t].size(); ++i)
      difference = std::max(difference, std::abs(sums_[t][i] - other.sums_[t][i]));
    return difference;
  }
  void add_man(int group, double born, double delta) {
    const size_t i = cell(group, born);
    men_[i] += delta;
    const int k = static_cast<int>(i % bins_);
    if (k < included_) {
      for (size_t t = 0; t < tables_.size(); ++t) spread(t, group, k, delta);
      if (++updates_ >= kPoolSumRefresh) refresh();
    }
  }
  void add_woman(int group, double born, double delta) { women_[cell(group, born)] += delta; }
  size_t bin_index(int group, double born) const { return cell(group, born); }
  size_t cells() const { return intensities_.size(); }
  // Cell-accounted fertility (D101) supplies channel totals directly before a
  // fertility-driven solve; per-woman contributions to that channel are zero.
  void assign_channel(int channel, const std::vector<double>& values) {
    for (size_t i = 0; i < intensities_.size(); ++i) intensities_[i][channel] = values[i];
  }
  void replace_role_intensity(int group, double born, const std::array<double, 3>& old,
                         const std::array<double, 3>& value) {
    const size_t i = cell(group, born);
    for (int c = 0; c < 3; ++c) {
      positive_contributors_[i][c] += (value[c] > 0) - (old[c] > 0);
      intensities_[i][c] += value[c] - old[c];
      // Exact zero when the last positive contribution leaves: subtraction
      // roundoff must not create phantom searchers in zero-risk cells.
      if (positive_contributors_[i][c] == 0) intensities_[i][c] = 0;
    }
  }
  void replace_intensity(int group, double born, const std::array<double, 2>& old,
                         const std::array<double, 2>& value) {
    replace_role_intensity(group, born, {old[0], old[1], 0}, {value[0], value[1], 0});
  }
  NumericMatrix intensity_snapshot(bool all_roles = false) const {
    const int channels = all_roles ? 3 : 2;
    NumericMatrix out(intensities_.size(), channels);
    for (size_t i = 0; i < intensities_.size(); ++i)
      for (int c = 0; c < channels; ++c) out(i, c) = intensities_[i][c];
    return out;
  }

  // Expected male flow shares follow from each woman-bin's intensity and
  // her group-level choice shares b[g] S[w,g] / sum_h b[h] S[w,h]. Find log
  // weights so those shares meet the mirrored target (Newton on a convex
  // objective, warm-started). Empty or uninformative states return b = 1.
  WeightSolution solve(double now, int channel, int table_id,
                       const std::vector<double>& delta, std::vector<double>& warm) {
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
    struct Row { int group, bin; double lambda; };
    std::vector<Row> rows;
    for (int g = 0; g < G; ++g) {
      for (int k = 0; k < bins_; ++k) {
        const double age = now - centre(k);
        if (age < cutoff_) continue;
        const size_t i = static_cast<size_t>(g) * bins_ + k;
        if (men_[i] > 0) { men[g].push_back({k, men_[i]}); men_count[g] += men_[i]; }
        if (women_[i] > 0) {
          women_count[g] += women_[i];
          const double lambda = std::max(0.0, intensities_[i][channel]);
          if (lambda > 0) { rows.push_back({g, k, lambda}); intensity[g] += lambda; }
        }
      }
    }
    double total_women = 0, total_intensity = 0;
    for (int g = 0; g < G; ++g) {
      total_women += women_count[g];
      total_intensity += intensity[g];
    }
    if (rows.empty() || !(total_women > 0) || !(total_intensity > 0)) { warm = s.log_weights; return s; }
    // Mirror: men of group g get the per-capita intensity of group-g women
    // (population-wide rate if none are present), then apply target deltas.
    const double overall = total_intensity / total_women;
    std::vector<double> mirror(G, 0);
    double mirror_sum = 0;
    for (int g = 0; g < G; ++g) {
      mirror[g] = men_count[g] * (women_count[g] > 0 ? intensity[g] / women_count[g] : overall);
      mirror_sum += mirror[g];
    }
    if (!(mirror_sum > 0)) { warm = s.log_weights; return s; }
    // Normalize log mirrored flows plus delta, preserving zero baseline targets.
    double largest = -std::numeric_limits<double>::infinity();
    for (int g = 0; g < G; ++g)
      if (mirror[g] > 0) largest = std::max(largest, std::log(mirror[g]) + delta[g]);
    double adjusted_sum = 0;
    for (int g = 0; g < G; ++g) {
      s.target[g] = mirror[g] > 0 ? std::exp(std::log(mirror[g]) + delta[g] - largest) : 0;
      adjusted_sum += s.target[g];
    }
    for (int g = 0; g < G; ++g) s.target[g] /= adjusted_sum;
    // Pool sums S[row, male group] come from the maintained sums. Entries below
    // the floor are recomputed exactly; a row that underflows is rescaled in
    // log space.
    const int R = static_cast<int>(rows.size());
    std::vector<double> S(static_cast<size_t>(R) * G, 0);
    for (int r = 0; r < R; ++r) {
      double total = 0;
      for (int m = 0; m < G; ++m) {
        double sum = sums[(static_cast<size_t>(rows[r].group) * G + m) * bins_ + rows[r].bin];
        if (!(sum >= kPoolSumFloor)) {
          sum = 0;
          for (const auto& bin : men[m])
            sum += bin.second * table.weight[table.index(rows[r].group, m, rows[r].bin - bin.first + bins_ - 1)];
        }
        S[static_cast<size_t>(r) * G + m] = sum;
        total += sum;
      }
      if (total > 0) continue;
      double maximum = R_NegInf;
      for (int m = 0; m < G; ++m)
        for (const auto& bin : men[m])
          maximum = std::max(maximum, table.utility[table.index(rows[r].group, m, rows[r].bin - bin.first + bins_ - 1)]);
      if (maximum == R_NegInf) continue;
      for (int m = 0; m < G; ++m) {
        double sum = 0;
        for (const auto& bin : men[m]) {
          const double u = table.utility[table.index(rows[r].group, m, rows[r].bin - bin.first + bins_ - 1)];
          if (u != R_NegInf) sum += bin.second * std::exp(u - maximum);
        }
        S[static_cast<size_t>(r) * G + m] = sum;
      }
    }
    std::vector<double> demand(G, 0);
    for (int r = 0; r < R; ++r)
      for (int m = 0; m < G; ++m) demand[m] += rows[r].lambda * S[static_cast<size_t>(r) * G + m];
    std::vector<int> free_groups;
    for (int g = 0; g < G; ++g) {
      if (!(s.target[g] > 0)) s.log_weights[g] = kZeroTargetLogWeight;
      else if (!(demand[g] > 0)) { s.log_weights[g] = kMaxLogWeight; s.infeasible = true; }
      else free_groups.push_back(g);
    }
    const int n = static_cast<int>(free_groups.size());
    if (n > 1) {
      double target_sum = 0;
      for (int g : free_groups) target_sum += s.target[g];
      std::vector<double> tau(n), x(n, 0);
      for (int j = 0; j < n; ++j) tau[j] = s.target[free_groups[j]] / target_sum;
      if (static_cast<int>(warm.size()) == G) {
        for (int j = 0; j < n; ++j) {
          const double start = warm[free_groups[j]] - warm[free_groups[0]];
          if (R_finite(start) && std::abs(start) < kMaxLogWeight) x[j] = start;
        }
      }
      x[0] = 0;
      std::vector<double> share(n), phi, hess, weights(n);
      auto evaluate = [&](const std::vector<double>& point, std::vector<double>& out_phi,
                          std::vector<double>& out_hess) {
        out_phi.assign(n, 0);
        out_hess.assign(static_cast<size_t>(n) * n, 0);
        double objective = 0, lambda_total = 0;
        for (int j = 0; j < n; ++j) weights[j] = std::exp(point[j]);
        for (int r = 0; r < R; ++r) {
          double denominator = 0;
          for (int j = 0; j < n; ++j) {
            share[j] = weights[j] * S[static_cast<size_t>(r) * G + free_groups[j]];
            denominator += share[j];
          }
          if (!(denominator > 0)) continue;
          const double lambda = rows[r].lambda;
          lambda_total += lambda;
          objective += lambda * std::log(denominator);
          for (int j = 0; j < n; ++j) {
            const double p = share[j] / denominator;
            out_phi[j] += lambda * p;
            for (int k = 0; k < n; ++k)
              out_hess[static_cast<size_t>(j) * n + k] +=
                lambda * ((j == k ? p : 0) - p * share[k] / denominator);
          }
        }
        for (double& value : out_phi) value /= lambda_total;
        for (double& value : out_hess) value /= lambda_total;
        objective /= lambda_total;
        for (int j = 0; j < n; ++j) objective -= tau[j] * point[j];
        return objective;
      };
      bool converged = false;
      std::vector<double> trial(n), trial_phi, trial_hess, step;
      for (int iteration = 0; iteration < 100; ++iteration) {
        const double objective = evaluate(x, phi, hess);
        double error = 0;
        for (int j = 0; j < n; ++j) error = std::max(error, std::abs(phi[j] - tau[j]));
        if (error < 1e-10) { converged = true; break; }
        // Newton step over the free log weights (the first is the reference).
        const int m = n - 1;
        std::vector<double> a(static_cast<size_t>(m) * m), b(m);
        for (int j = 0; j < m; ++j) {
          b[j] = -(phi[j + 1] - tau[j + 1]);
          for (int k = 0; k < m; ++k)
            a[static_cast<size_t>(j) * m + k] = hess[static_cast<size_t>(j + 1) * n + k + 1];
          a[static_cast<size_t>(j) * m + j] += 1e-14;
        }
        if (!solve_linear(a, b, m, step)) break;
        // Accept a step that lowers the convex objective or the target mismatch;
        // near the solution, roundoff can hide the objective's decrease.
        double length = 1;
        bool moved = false;
        while (length > 1e-12) {
          trial[0] = 0;
          for (int j = 0; j < m; ++j)
            trial[j + 1] = std::max(-kMaxLogWeight, std::min(kMaxLogWeight, x[j + 1] + length * step[j]));
          const double trial_objective = evaluate(trial, trial_phi, trial_hess);
          double trial_error = 0;
          for (int j = 0; j < n; ++j) trial_error = std::max(trial_error, std::abs(trial_phi[j] - tau[j]));
          if (trial_objective < objective || trial_error < error) { moved = true; break; }
          length /= 2;
        }
        if (!moved) break;
        x = trial;
      }
      if (!converged) s.infeasible = true;
      for (int j = 0; j < n; ++j) s.log_weights[free_groups[j]] = x[j];
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
        const double pool = S[static_cast<size_t>(r) * G + m];
        if (pool > 0) top = std::max(top, s.log_weights[m] + std::log(pool));
      }
      if (top == R_NegInf) continue;
      double denominator = 0;
      for (int m = 0; m < G; ++m) {
        const double pool = S[static_cast<size_t>(r) * G + m];
        share[m] = pool > 0 ? std::exp(s.log_weights[m] + std::log(pool) - top) : 0;
        denominator += share[m];
      }
      flow_total += rows[r].lambda;
      for (int m = 0; m < G; ++m) s.flow[m] += rows[r].lambda * share[m] / denominator;
    }
    if (flow_total > 0) {
      double miss = 0;
      for (int m = 0; m < G; ++m) {
        s.flow[m] /= flow_total;
        miss = std::max(miss, std::abs(s.flow[m] - s.target[m]));
      }
      if (miss > 1e-8) s.infeasible = true;
    }
    warm = s.log_weights;
    return s;
  }
};

// Actual current hazard and next absolute age boundary, including schedule gaps.
// Absolute comparisons make a queued boundary use the right-hand band even
// when subtracting birth time would round the attained age just below its edge.
std::pair<double, double> current_hazard(const std::vector<Band>* bands,
                                         double born, double now) {
  double value = 0, next = R_PosInf;
  if (bands) for (const Band& band : *bands) {
    const double start = born + band.start, end = born + band.end;
    if (start <= now && now < end) value = band.rate;
    if (start > now) next = std::min(next, start);
    if (end > now) next = std::min(next, end);
  }
  return {value, next};
}

// Ordered bookkeeping boundaries with at most one pending time per person.
// Shared by the calibration tracker and fertility cells: entries are rate
// changes (age bands, gaps, spacing release), not demographic events, and
// consume no randomness. Ties are ordered by person index.
class BoundaryQueue {
  // Binary min-heap of (time, person) with lazy deletion: an entry is current
  // only while it matches next_[person]. Pops follow (time, person) order.
  std::vector<std::pair<double, int>> heap_;
  std::vector<double> next_;
  size_t live_ = 0;
  struct Later {
    bool operator()(const std::pair<double, int>& a, const std::pair<double, int>& b) const {
      return a > b;
    }
  };
  void discard_stale() {
    while (!heap_.empty() && heap_.front().first != next_[heap_.front().second]) {
      std::pop_heap(heap_.begin(), heap_.end(), Later());
      heap_.pop_back();
    }
  }
 public:
  void set(int person, double time) {
    if (person >= static_cast<int>(next_.size())) next_.resize(person + 1, R_PosInf);
    if (next_[person] == time) return;
    if (R_finite(next_[person])) --live_;
    next_[person] = time;
    if (R_finite(time)) {
      ++live_;
      heap_.push_back({time, person});
      std::push_heap(heap_.begin(), heap_.end(), Later());
    }
    // Bound stale entries: rebuild from current times when they dominate.
    if (heap_.size() > 64 && heap_.size() > 4 * live_) {
      heap_.clear();
      for (size_t p = 0; p < next_.size(); ++p)
        if (R_finite(next_[p])) heap_.push_back({next_[p], static_cast<int>(p)});
      std::make_heap(heap_.begin(), heap_.end(), Later());
    }
  }
  void clear(int person) {
    if (person < static_cast<int>(next_.size())) set(person, R_PosInf);
  }
  double top() {
    discard_stale();
    return heap_.empty() ? R_PosInf : heap_.front().first;
  }
  // Remove and return the earliest entry; the caller refreshes that person.
  std::pair<double, int> pop() {
    discard_stale();
    const auto entry = heap_.front();
    std::pop_heap(heap_.begin(), heap_.end(), Later());
    heap_.pop_back();
    next_[entry.second] = R_PosInf;
    --live_;
    return entry;
  }
  size_t size() const { return live_; }
};

// One pending boundary per active woman, replaced eagerly on state changes.
// Separate channels share one queue and one state-update interface; no Cartesian
// table of previous status, parity, region, etc. is needed by calibration.
class FemaleIntensityTracker {
  struct Woman {
    int group = 0;
    double born = 0, search_origin = 0, ready = 0;
    bool active = false;
    const std::vector<Band>* search = nullptr;
    const std::vector<Band>* fertility = nullptr;
    const std::vector<Band>* probability = nullptr;
    std::array<double, 3> value = {0, 0, 0};
  };
  WeightCalibrator& calibration_;
  std::vector<Woman> women_;
  BoundaryQueue queue_;
  unsigned long long updates_ = 0;
  void refresh(int person, double now) {
    Woman& w = women_[person];
    const auto search = current_hazard(w.search, w.search_origin, now);
    const auto fertility = now < w.ready ? std::make_pair(0.0, w.ready) :
      current_hazard(w.fertility, w.born, now);
    const auto probability = current_hazard(w.probability, w.born, now);
    const double p = probability.first;
    const std::array<double, 3> value = {search.first * (1 - p), fertility.first, search.first * p};
    calibration_.replace_role_intensity(w.group, w.born, w.value, value);
    w.value = value;
    queue_.set(person, std::min(std::min(search.second, fertility.second), probability.second));
  }
 public:
  explicit FemaleIntensityTracker(WeightCalibrator& calibration) : calibration_(calibration) {}
  void activate(int person, int group, double born, const std::vector<Band>* search,
                const std::vector<Band>* fertility, double ready, double now,
                double search_origin = NA_REAL, const std::vector<Band>* probability = nullptr) {
    if (person >= static_cast<int>(women_.size())) women_.resize(person + 1);
    remove(person);
    Woman& w = women_[person];
    w.active = true; w.group = group; w.born = born;
    w.search_origin = NumericVector::is_na(search_origin) ? born : search_origin;
    w.search = search; w.fertility = fertility; w.ready = ready; w.probability = probability;
    refresh(person, now);
  }
  void remove(int person) {
    if (person >= static_cast<int>(women_.size()) || !women_[person].active) return;
    Woman& w = women_[person];
    queue_.clear(person);
    calibration_.replace_role_intensity(w.group, w.born, w.value, {0, 0, 0});
    w.active = false; w.value = {0, 0, 0};
  }
  void set_ready(int person, double ready, double now) {
    if (person >= static_cast<int>(women_.size()) || !women_[person].active) return;
    women_[person].ready = ready;
    refresh(person, now);
  }
  void advance(double now) {
    while (queue_.top() <= now) {
      const auto entry = queue_.pop();
      refresh(entry.second, entry.first);
      if (++updates_ % 1024 == 0) checkUserInterrupt();
    }
  }
  size_t pending() const { return queue_.size(); }
};

// ---- Cell-accounted fertility under birth spacing (D101) -------------------
// A cell is one age sub-band (at most one year, split at whole ages and input
// band edges) of one effective fertility schedule (group, union type, and
// actual parity when schedules are parity-specific; pooled bands are split by
// actual parity). Each cell generates births at m_c * N_c, counting
// all living members including those inside the birth interval, and assigns
// each birth to an eligible member, weighted by persistent hetfert when
// enabled (selection probability, not a rate multiplier). Without hetfert,
// selection is uniform and no weight tree is allocated. A cell with members but no eligible member
// cannot produce births: its intensity is dropped and recorded as unrealized.
// Clocks keep an Exp(1) threshold consumed at the current intensity, so count
// changes recompute firing times without drawing random numbers.
// Upper-open fertility bands stay one cell when birth_interval * rate is below
// this share (a numerical tolerance, not a demographic parameter).
constexpr double kOpenBandSplitTolerance = 0.001;

// Positive eligible weights in swap-remove slots. Recompute sums from children
// rather than subtracting from a running total, so removing a dominant weight
// restores the small weights accurately. Capacity grows geometrically.
class FertilityWeightTree {
  size_t size_ = 0, capacity_ = 0;
  std::vector<double> tree_;
  void set(size_t slot, double value) {
    size_t node = capacity_ + slot;
    tree_[node] = value;
    while (node > 1) {
      node /= 2;
      tree_[node] = tree_[2 * node] + tree_[2 * node + 1];
      if (!R_finite(tree_[node])) stop("Eligible hetfert weight sum is not representable.");
    }
  }
 public:
  double sum() const { return tree_.empty() ? 0 : tree_[1]; }
  void add(double value) {
    if (size_ == capacity_) {
      const size_t next = capacity_ ? 2 * capacity_ : 1;
      std::vector<double> grown(2 * next, 0);
      for (size_t i = 0; i < size_; ++i) grown[next + i] = tree_[capacity_ + i];
      for (size_t i = next - 1; i > 0; --i) grown[i] = grown[2 * i] + grown[2 * i + 1];
      capacity_ = next;
      tree_.swap(grown);
    }
    set(size_++, value);
  }
  void remove(size_t slot) {
    const size_t last = size_ - 1;
    const double moved = tree_[capacity_ + last];
    set(slot, 0);
    if (slot != last) { set(last, 0); set(slot, moved); }
    --size_;
  }
  size_t choose(double uniform) const {
    double value = std::min(uniform * sum(), std::nextafter(sum(), 0.0));
    size_t node = 1;
    while (node < capacity_) {
      const double left = tree_[2 * node];
      if (value < left) node *= 2;
      else { value -= left; node = 2 * node + 1; }
    }
    return std::min(node - capacity_, size_ - 1);
  }
};

// Log-space form of the accepted latent-normal inheritance; valid also at
// sigma=0 for supplied nonunit weights. Endpoints avoid needless RNG draws.
double inherit_fertility_weight(double maternal, double variance, double alpha) {
  if (alpha == 1) return maternal;
  const double mu = -variance / 2;
  double log_weight = mu + alpha * (std::log(maternal) - mu);
  if (variance > 0) log_weight += std::sqrt(variance) * std::sqrt(1 - alpha * alpha) * R::rnorm(0, 1);
  const double weight = std::exp(log_weight);
  if (!(weight > 0) || !R_finite(weight)) stop("Generated hetfert weight is not finite and strictly positive.");
  return weight;
}

// Compensated bin sums retain small weights when a dominant woman leaves.
struct FertilityWeightSum {
  double sum = 0, correction = 0;
  void add(double value) {
    const double next = sum + value;
    correction += std::abs(sum) >= std::abs(value) ? (sum - next) + value : (value - next) + sum;
    sum = next;
  }
  double value() const { return sum + correction; }
};

class FertilityCells {
 public:
  struct State {
    int stratum = -1;          // -1: not a living woman
    double born = 0, ready = 0, weight = 1;
    bool unpartnered = false;
  };
  // An input schedule for one group/type and input parity band.
  struct Source {
    int group = 0, type = 0;
    double parity_start = 0, parity_end = R_PosInf;
    std::vector<Band> bands;
  };
  // Cells of one group/type and actual parity (or all parities when schedules
  // are not parity-specific), with age bands split to at most one year.
  struct Stratum {
    int group = 0, type = 0;
    double parity_start = 0, parity_end = R_PosInf;
    std::vector<Band> bands;
    int first_cell = 0;
  };

 private:
  struct Cell {
    int stratum = 0;
    double rate = 0, age_start = 0, age_end = 0;
    int members = 0;
    std::vector<int> eligible;                 // swap-remove list of women
    FertilityWeightTree weights;
    double lambda = 0, threshold = 0, last = 0, due = R_PosInf;
    bool has_threshold = false, queued = false;
    std::map<size_t, std::pair<int, FertilityWeightSum>> counted; // count and weight per unpartnered calibration bin
  };
  struct Member {
    int cell = -1, slot = -1;
    bool eligible = false, counted = false;
    size_t bin = 0;
    double weight = 1;
  };
  std::vector<Source> sources_;
  std::vector<Stratum> strata_;
  std::vector<Cell> cells_;
  std::vector<Member> members_;
  // lookup_[group * types + type][actual parity] -> stratum (-1 until created);
  // without parity-specific schedules, index 0 serves every parity.
  int types_;
  bool split_parity_, weighted_;
  double max_age_, interval_;
  std::vector<std::vector<int>> lookup_;
  BoundaryQueue boundaries_;
  // Cell birth clocks: an indexed binary min-heap of armed cells keyed by
  // (due, cell), updated in place without allocation. The main loop compares
  // its head with the person event queue.
  std::vector<int> heap_, position_;
  bool earlier(int a, int b) const {
    const double da = cells_[a].due, db = cells_[b].due;
    return da < db || (da == db && a < b);
  }
  void swap_nodes(size_t i, size_t j) {
    std::swap(heap_[i], heap_[j]);
    position_[heap_[i]] = static_cast<int>(i);
    position_[heap_[j]] = static_cast<int>(j);
  }
  void sift_up(size_t i) {
    while (i > 0) {
      const size_t parent = (i - 1) / 2;
      if (!earlier(heap_[i], heap_[parent])) break;
      swap_nodes(i, parent);
      i = parent;
    }
  }
  void sift_down(size_t i) {
    for (;;) {
      size_t best = i;
      const size_t left = 2 * i + 1, right = left + 1;
      if (left < heap_.size() && earlier(heap_[left], heap_[best])) best = left;
      if (right < heap_.size() && earlier(heap_[right], heap_[best])) best = right;
      if (best == i) break;
      swap_nodes(i, best);
      i = best;
    }
  }
  void unschedule(int c) {
    const size_t i = static_cast<size_t>(position_[c]);
    const size_t last = heap_.size() - 1;
    if (i != last) swap_nodes(i, last);
    heap_.pop_back();
    position_[c] = -1;
    if (i < heap_.size()) { sift_up(i); sift_down(i); }
  }
  void schedule_clock(int c) {
    if (static_cast<size_t>(c) >= position_.size()) position_.resize(c + 1, -1);
    position_[c] = static_cast<int>(heap_.size());
    heap_.push_back(c);
    sift_up(heap_.size() - 1);
  }
  std::function<State(int)> state_;
  std::function<double()> draw_threshold_;
  const WeightCalibrator* calibrator_;
  double start_, end_;
  // Unrealized expected births and member-years without eligible women, by
  // one-year period from the simulation start and cell.
  std::map<std::pair<long long, int>, std::array<double, 2>> unrealized_;
  unsigned long long updates_ = 0;

  void record_unrealized(int c, double from, double to, double lambda) {
    while (from < to) {
      long long k = static_cast<long long>(std::floor(from - start_));
      if (start_ + static_cast<double>(k + 1) <= from) ++k; // roundoff at a boundary
      const double stop = std::min(to, start_ + static_cast<double>(k + 1));
      auto& value = unrealized_[{k, c}];
      value[0] += lambda * (stop - from);
      value[1] += stop - from;
      from = stop;
    }
  }
  // Account for the interval since the last update under the state that held
  // during it. Must run before membership or eligibility changes.
  void settle(int c, double now) {
    Cell& x = cells_[c];
    if (x.lambda > 0 && now > x.last) {
      if (x.eligible.empty()) record_unrealized(c, x.last, now, x.lambda);
      else if (x.has_threshold) x.threshold = std::max(0.0, x.threshold - x.lambda * (now - x.last));
    }
    x.last = now;
  }
  // Recompute the intensity and firing time from current counts (after settle).
  void touch(int c, double now) {
    Cell& x = cells_[c];
    x.lambda = x.rate * x.members;
    rearm(c, now);
  }
  void rearm(int c, double now) {
    Cell& x = cells_[c];
    if (x.queued) { unschedule(c); x.queued = false; }
    if (!(x.lambda > 0) || x.eligible.empty()) return;
    if (!x.has_threshold) { x.threshold = draw_threshold_(); x.has_threshold = true; }
    const double due = now + x.threshold / x.lambda;
    if (due < end_) {
      x.due = due;
      schedule_clock(c);
      x.queued = true;
    }
  }
  void leave(int woman, int& changed) {
    Member& m = members_[woman];
    if (m.cell < 0) return;
    Cell& x = cells_[m.cell];
    --x.members;
    if (m.eligible) {
      const int last = x.eligible.back();
      x.eligible[m.slot] = last;
      members_[last].slot = m.slot;
      x.eligible.pop_back();
      if (weighted_) x.weights.remove(m.slot);
    }
    if (m.counted) {
      auto& bin = x.counted[m.bin];
      if (weighted_) bin.second.add(-m.weight);
      if (--bin.first == 0) x.counted.erase(m.bin);
    }
    changed = m.cell;
    m = Member();
  }

 public:
  // Split input age bands at whole ages and band edges. An upper-open band is
  // split up to the oldest attainable age in the run, then left open, unless
  // interval * rate is below kOpenBandSplitTolerance: at most that share of its
  // members can be in spacing, so splitting could not change recovery
  // materially and would only add yearly membership updates.
  static std::vector<Band> split_ages(const std::vector<Band>& bands, double max_age,
                                      double interval) {
    std::vector<Band> out;
    for (const Band& band : bands) {
      // Zero-rate bands generate no births and need no members: treat as gaps.
      if (!(band.rate > 0)) continue;
      double lo = band.start;
      if (!std::isfinite(band.end) && interval * band.rate < kOpenBandSplitTolerance) {
        out.push_back(band);
        continue;
      }
      const double cap = std::isfinite(band.end) ? band.end : std::max(lo, std::ceil(max_age));
      while (lo < cap) {
        const double next = std::min(cap, std::floor(lo) + 1);
        out.push_back({lo, next, band.rate});
        lo = next;
      }
      if (!std::isfinite(band.end)) out.push_back({lo, band.end, band.rate});
    }
    return out;
  }
  int add_stratum(const Source& source, double parity_start, double parity_end) {
    Stratum st;
    st.group = source.group; st.type = source.type;
    st.parity_start = parity_start; st.parity_end = parity_end;
    st.bands = split_ages(source.bands, max_age_, interval_);
    st.first_cell = static_cast<int>(cells_.size());
    const int k = static_cast<int>(strata_.size());
    for (const Band& band : st.bands) {
      Cell x;
      x.stratum = k;
      x.rate = band.rate; x.age_start = band.start; x.age_end = band.end;
      x.last = start_;
      cells_.push_back(x);
    }
    strata_.push_back(st);
    return k;
  }

 public:
  FertilityCells(std::vector<Source> sources, int groups, int types, bool split_parity,
                 double max_age, double interval, bool weighted,
                 std::function<State(int)> state, std::function<double()> draw_threshold,
                 const WeightCalibrator* calibrator, double start, double end)
      : sources_(std::move(sources)), types_(types), split_parity_(split_parity), weighted_(weighted),
        max_age_(max_age), interval_(interval), lookup_(static_cast<size_t>(groups) * types),
        state_(std::move(state)), draw_threshold_(std::move(draw_threshold)),
        calibrator_(calibrator), start_(start), end_(end) {
    if (!split_parity_)
      for (const Source& source : sources_)
        lookup_[static_cast<size_t>(source.group) * types_ + source.type] =
          {add_stratum(source, source.parity_start, source.parity_end)};
  }
  // Input parity bands partition [0, Inf) within each group/type (validated in
  // R). Actual-parity strata are created when a woman first reaches a parity.
  int stratum(int group, int type, int parity) {
    auto& list = lookup_[static_cast<size_t>(group) * types_ + type];
    if (!split_parity_) return list.empty() ? -1 : list[0];
    if (parity < 0) return -1;
    if (parity < static_cast<int>(list.size()) && list[parity] >= 0) return list[parity];
    for (const Source& source : sources_) {
      if (source.group != group || source.type != type ||
          !(source.parity_start <= parity && parity < source.parity_end)) continue;
      if (parity >= static_cast<int>(list.size())) list.resize(parity + 1, -1);
      list[parity] = add_stratum(source, parity, parity + 1.0);
      return list[parity];
    }
    return -1;
  }
  // Recompute a woman's cell, eligibility and calibration bin from engine state.
  void refresh(int woman, double now) {
    if (woman >= static_cast<int>(members_.size())) members_.resize(woman + 1);
    const State s = state_(woman);
    int cell = -1;
    double next = R_PosInf;
    if (s.stratum >= 0) {
      // Bands are sorted and disjoint: find the last band starting at or
      // before now, using absolute times as in current_hazard().
      const Stratum& st = strata_[s.stratum];
      const auto& bands = st.bands;
      size_t lo = 0, hi = bands.size();
      while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (s.born + bands[mid].start <= now) lo = mid + 1; else hi = mid;
      }
      if (lo > 0 && now < s.born + bands[lo - 1].end) {
        cell = st.first_cell + static_cast<int>(lo - 1);
        next = s.born + bands[lo - 1].end;
      } else if (lo < bands.size()) {
        next = s.born + bands[lo].start;
      }
    }
    const bool eligible = cell >= 0 && now >= s.ready;
    if (cell >= 0 && !eligible) next = std::min(next, s.ready);
    const bool counted = calibrator_ && eligible && s.unpartnered;
    const size_t bin = counted ? calibrator_->bin_index(strata_[s.stratum].group, s.born) : 0;
    Member& m = members_[woman];
    if (m.cell != cell || m.eligible != eligible || m.counted != counted || (counted && m.bin != bin)) {
      // Settle elapsed intensity under the old state before changing it.
      if (m.cell >= 0) settle(m.cell, now);
      if (cell >= 0 && cell != m.cell) settle(cell, now);
      int old_cell = -1;
      leave(woman, old_cell);
      if (cell >= 0) {
        Cell& x = cells_[cell];
        ++x.members;
        m.cell = cell; m.eligible = eligible; m.counted = counted; m.bin = bin; m.weight = s.weight;
        if (eligible) {
          m.slot = static_cast<int>(x.eligible.size()); x.eligible.push_back(woman);
          if (weighted_) x.weights.add(m.weight);
        }
        if (counted) { ++x.counted[bin].first; if (weighted_) x.counted[bin].second.add(m.weight); }
      }
      if (old_cell >= 0) touch(old_cell, now);
      if (cell >= 0 && cell != old_cell) touch(cell, now);
    }
    boundaries_.set(woman, next);
  }
  double next_boundary() { return boundaries_.top(); }
  void advance_boundary() {
    const auto entry = boundaries_.pop();
    refresh(entry.second, entry.first);
    if (++updates_ % 1024 == 0) checkUserInterrupt();
  }
  double next_clock() const { return heap_.empty() ? R_PosInf : cells_[heap_.front()].due; }
  // The earliest cell clock fires: spend its threshold, pick a mother, re-arm.
  int fire_next(double now) {
    const int c = heap_.front();
    unschedule(c);
    Cell& x = cells_[c];
    x.queued = false;
    x.has_threshold = false;
    x.last = now;
    if (x.eligible.empty()) stop("A fertility cell fired without eligible women.");
    const double n = static_cast<double>(x.eligible.size());
    const double uniform = R::runif(0, 1);
    const size_t pick = weighted_ ? x.weights.choose(uniform) :
      std::min(x.eligible.size() - 1, static_cast<size_t>(std::floor(uniform * n)));
    const int woman = x.eligible[pick];
    rearm(c, now);
    return woman;
  }
  // Calibration fertility totals: each unpartnered eligible member contributes
  // her equivalent weighted hazard in her calibration bin.
  std::vector<double> fill_calibration(WeightCalibrator& calibrator) const {
    std::vector<double> values(calibrator.cells(), 0.0);
    for (const Cell& x : cells_) {
      if (!(x.lambda > 0) || x.eligible.empty() || x.counted.empty()) continue;
      if (weighted_) {
        const double total = x.weights.sum();
        for (const auto& bin : x.counted)
          values[bin.first] += x.lambda * (bin.second.second.value() / total);
      } else {
        const double hazard = x.lambda / static_cast<double>(x.eligible.size());
        for (const auto& bin : x.counted) values[bin.first] += hazard * bin.second.first;
      }
    }
    calibrator.assign_channel(1, values);
    return values;
  }
  // A woman's equivalent individual hazard (rule 4), zero inside spacing.
  double hazard(int woman) const {
    if (woman >= static_cast<int>(members_.size())) return 0;
    const Member& m = members_[woman];
    if (m.cell < 0 || !m.eligible) return 0;
    const Cell& x = cells_[m.cell];
    return weighted_ ? x.lambda * (m.weight / x.weights.sum()) :
      x.lambda / static_cast<double>(x.eligible.size());
  }
  // Close unrealized-intensity intervals at the end of the run.
  void finish(double end) {
    for (size_t c = 0; c < cells_.size(); ++c) {
      Cell& x = cells_[c];
      if (x.lambda > 0 && x.eligible.empty() && end > x.last)
        record_unrealized(static_cast<int>(c), x.last, end, x.lambda);
      x.last = end;
    }
  }
  // Internal consistency check against an explicit scan of member records.
  bool consistent() const {
    std::vector<int> members(cells_.size(), 0), eligible(cells_.size(), 0);
    std::vector<double> weights(cells_.size(), 0);
    for (size_t w = 0; w < members_.size(); ++w) {
      const Member& m = members_[w];
      if (m.cell < 0) continue;
      ++members[m.cell];
      if (m.eligible) {
        ++eligible[m.cell];
        weights[m.cell] += m.weight;
        if (cells_[m.cell].eligible[m.slot] != static_cast<int>(w)) return false;
      }
    }
    for (size_t c = 0; c < cells_.size(); ++c)
      if (members[c] != cells_[c].members || eligible[c] != static_cast<int>(cells_[c].eligible.size()))
        return false;
    if (weighted_) for (size_t c = 0; c < cells_.size(); ++c)
      if (std::abs(weights[c] - cells_[c].weights.sum()) > 1e-9 * std::max(weights[c], 1.0)) return false;
    return true;
  }
  DataFrame unrealized(const std::vector<std::string>& type_labels, double cutoff = R_PosInf) const {
    const size_t n = unrealized_.size();
    NumericVector period_start(n), period_end(n), parity_start(n), parity_end(n),
      age_start(n), age_end(n), births(n), years(n);
    IntegerVector group(n);
    CharacterVector type(n);
    size_t r = 0;
    for (const auto& entry : unrealized_) {
      const Cell& x = cells_[entry.first.second];
      const Stratum& st = strata_[x.stratum];
      period_start[r] = start_ + static_cast<double>(entry.first.first);
      period_end[r] = std::min(std::min(end_, cutoff), period_start[r] + 1);
      group[r] = st.group;
      type[r] = type_labels[st.type];
      parity_start[r] = st.parity_start; parity_end[r] = st.parity_end;
      age_start[r] = x.age_start; age_end[r] = x.age_end;
      births[r] = entry.second[0]; years[r] = entry.second[1];
      ++r;
    }
    return DataFrame::create(_["period_start"] = period_start, _["period_end"] = period_end,
      _["group"] = group, _["union_type"] = type, _["parity_start"] = parity_start,
      _["parity_end"] = parity_end, _["age_start"] = age_start, _["age_end"] = age_end,
      _["unrealized_births"] = births, _["years_without_eligible"] = years);
  }
};

// Input schedules for cells from the effective (group-multiplied) fertility table.
std::vector<FertilityCells::Source> fertility_strata(DataFrame rates, int groups,
    const std::vector<std::string>& type_labels) {
  std::vector<FertilityCells::Source> out;
  const bool parity = rates.containsElementNamed("parity_start");
  const bool grouped = rates.containsElementNamed(".hammel_group");
  for (int g = 0; g < groups; ++g) {
    for (int t = 0; t < static_cast<int>(type_labels.size()); ++t) {
      const std::string label = type_labels[t];
      std::map<double, double> bands = {{0, R_PosInf}};
      if (parity) {
        bands.clear();
        NumericVector starts = rates["parity_start"], ends = rates["parity_end"];
        IntegerVector row_groups;
        if (grouped) row_groups = rates[".hammel_group"];
        CharacterVector row_types;
        if (!label.empty()) row_types = rates["union_type"];
        for (int i = 0; i < starts.size(); ++i) {
          if (grouped && row_groups[i] != g) continue;
          if (!label.empty() && as<std::string>(row_types[i]) != label) continue;
          bands[starts[i]] = ends[i];
        }
      }
      for (const auto& band : bands) {
        FertilityCells::Source st;
        st.group = g; st.type = t;
        st.parity_start = band.first; st.parity_end = band.second;
        st.bands = read_bands(rates, "", label, g, "", false, parity ? band.first : -1);
        out.push_back(st);
      }
    }
  }
  return out;
}

// Internal entry point: calibrate weights for a supplied unpartnered state.
// [[Rcpp::export]]
List inspect_weighted_calibration_cpp(NumericVector female_born, IntegerVector female_groups,
                                      LogicalVector female_fertile, NumericVector male_born,
                                      IntegerVector male_groups, DataFrame rates, List matching,
                                      int group_count, double now, double cutoff, double width,
                                      bool fertility,
                                      Nullable<IntegerVector> female_previous = R_NilValue,
                                      Nullable<NumericVector> female_previous_ends = R_NilValue) {
  if (female_groups.size() != female_born.size() || female_fertile.size() != female_born.size() ||
      male_groups.size() != male_born.size() || group_count < 1)
    stop("Invalid weighted-calibration inputs.");
  double lo = now, hi = now;
  for (double born : female_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  for (double born : male_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  const bool previous_specific = rates.containsElementNamed("previous_union") && !fertility;
  IntegerVector previous(female_born.size(), 0);
  if (female_previous.isNotNull()) previous = IntegerVector(female_previous);
  if (previous.size() != female_born.size()) stop("Invalid previous-union calibration inputs.");
  for (int p : previous) if (p < 0 || p > 2) stop("Invalid previous-union calibration category.");
  WeightCalibrator calibrator(lo, hi, group_count, width, cutoff);
  const LagTable table = build_lag_table(matching, group_count, calibrator.bins(), width);
  calibrator.add_table(&table);
  for (R_xlen_t i = 0; i < male_born.size(); ++i) calibrator.add_man(male_groups[i], male_born[i], 1);
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    calibrator.add_woman(female_groups[i], female_born[i], 1);
  }
  NumericVector previous_ends(female_born.size(), NA_REAL);
  if (female_previous_ends.isNotNull()) previous_ends = NumericVector(female_previous_ends);
  if (previous_ends.size() != female_born.size()) stop("Invalid previous ending times.");
  const int categories = previous_specific ? 3 : 1;
  std::vector<bool> duration_bands(group_count * categories, false);
  std::vector<std::vector<Band>> bands(group_count * categories);
  std::vector<const std::vector<Band>*> hazards(group_count * categories);
  const std::array<std::string, 3> labels = {"none", "separated", "widowed"};
  for (int g = 0; g < group_count; ++g) {
    for (int p = 0; p < categories; ++p) {
      const int i = g * categories + p;
      duration_bands[i] = previous_specific && search_uses_duration(rates, labels[p]);
      bands[i] = read_bands(rates, "", "", g, previous_specific ? labels[p] : "", duration_bands[i]);
      hazards[i] = &bands[i];
    }
  }
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    const int schedule = female_groups[i] * categories + (previous_specific ? previous[i] : 0);
    const double origin = duration_bands[schedule] ? previous_ends[i] : female_born[i];
    if (!R_finite(origin)) stop("Duration calibration requires a previous ending time.");
    const double value = current_hazard(hazards[schedule], origin, now).first;
    calibrator.replace_intensity(female_groups[i], female_born[i], {0, 0},
      {value, female_fertile[i] ? value : 0});
  }
  std::vector<double> warm;
  const WeightSolution s = calibrator.solve(now, fertility, 0,
    as<std::vector<double>>(matching["delta"]), warm);
  return List::create(_["log_weights"] = wrap(s.log_weights), _["target"] = wrap(s.target),
                      _["flow"] = wrap(s.flow), _["infeasible"] = s.infeasible);
}

// ---- Exact full-pool partner draws (D072) ----------------------------------
// Living, unpartnered men are kept in cells by group and birth-time bin. A draw
// picks a cell with probability proportional to its size times the largest
// weight any member could have for this woman, picks a member uniformly,
// accepts him with probability exp(U - cell maximum), and screens only him.
// This samples conditional-logit choice over the full eligible pool exactly,
// without a candidate sample.
class MalePool {
  int groups_, bins_;
  long long first_;
  double width_;
  std::vector<std::vector<int>> cells_;
  std::vector<int> slot_, cell_;
  size_t size_ = 0;

 public:
  MalePool(double min_birth, double max_birth, int groups, double width)
      : groups_(groups), width_(width) {
    if (!(width > 0) || !R_finite(min_birth) || !R_finite(max_birth) || max_birth < min_birth ||
        groups < 1)
      stop("Invalid male pool range.");
    first_ = static_cast<long long>(std::floor(min_birth / width));
    const long long last = static_cast<long long>(std::floor(max_birth / width));
    if (last - first_ + 1 > 50000000LL / groups) stop("Male pool birth-time range is too large.");
    bins_ = static_cast<int>(last - first_ + 1);
    cells_.assign(static_cast<size_t>(groups) * bins_, std::vector<int>());
  }
  int groups() const { return groups_; }
  int bins() const { return bins_; }
  long long first() const { return first_; }
  double width() const { return width_; }
  bool empty() const { return size_ == 0; }
  const std::vector<int>& members(int group, int bin) const {
    return cells_[static_cast<size_t>(group) * bins_ + bin];
  }
  bool contains(int person) const {
    return static_cast<size_t>(person) < slot_.size() && slot_[person] >= 0;
  }
  void add(int person, int group, double born) {
    if (static_cast<size_t>(person) >= slot_.size()) {
      slot_.resize(person + 1, -1);
      cell_.resize(person + 1, -1);
    }
    if (slot_[person] >= 0) return;
    const long long k = static_cast<long long>(std::floor(born / width_)) - first_;
    if (k < 0 || k >= bins_ || group < 0 || group >= groups_)
      stop("Person is outside the male pool range.");
    const int c = group * bins_ + static_cast<int>(k);
    slot_[person] = static_cast<int>(cells_[c].size());
    cell_[person] = c;
    cells_[c].push_back(person);
    ++size_;
  }
  void remove(int person) {
    if (!contains(person)) return;
    std::vector<int>& members = cells_[cell_[person]];
    const int moved = members.back();
    members[slot_[person]] = moved;
    slot_[moved] = slot_[person];
    members.pop_back();
    slot_[person] = -1;
    --size_;
  }
  // Last bin whose oldest possible member is at least `cutoff` years old; the
  // youngest such bin can hold younger men, who are rejected individually.
  int last_bin(double now, double cutoff) const {
    const double k = std::floor((now - cutoff) / width_) - static_cast<double>(first_);
    return static_cast<int>(std::min(static_cast<double>(bins_ - 1), std::max(-1.0, k)));
  }
};

// One matching model prepared for exact draws. The native linear/quadratic age
// model has a closed-form maximum over any age-gap interval. Custom formulas
// are bounded by their values on a fine gap grid (evaluated once); a drawn
// man's utility is always evaluated exactly, and a bound found too low is
// counted so simulation can report approximate draws.
struct DrawModel {
  List model;
  bool native = false;
  double linear = 0, quadratic = 0;
  NumericMatrix offsets;
  LogicalMatrix blocked;
  bool prevent_all = false;
  int groups = 0;
  double grid_start = 0, grid_step = 0;
  int grid_size = 0;
  std::vector<double> grid; // [(female * groups + male) * grid_size + j]

  double offset(int f, int m) const { return offsets.nrow() ? offsets(f, m) : 0.0; }
  bool allowed(int f, int m) const {
    return !(prevent_all && f != m) && !(blocked.nrow() && blocked(f, m));
  }
  double native_value(double gap) const { return linear * gap + quadratic * gap * gap; }
  // Exact utilities, including pair effects, for one female/male group pair.
  NumericVector utilities(int f, int m, const NumericVector& gaps) const {
    if (native) {
      NumericVector out(gaps.size());
      for (R_xlen_t j = 0; j < gaps.size(); ++j) {
        out[j] = native_value(gaps[j]) + offset(f, m);
        if (!R_finite(out[j])) stop("Matching ages and utilities must be finite.");
      }
      return out;
    }
    return matching_utilities(gaps, model, -1, IntegerVector(gaps.size(), m),
                              IntegerVector(gaps.size(), f));
  }
  // Upper bound on utility for male-minus-female age gaps in [lo, hi].
  double bound(int f, int m, double lo, double hi) const {
    if (native) {
      double best = std::max(native_value(lo), native_value(hi));
      if (quadratic < 0) {
        const double vertex = -linear / (2 * quadratic);
        if (vertex > lo && vertex < hi) best = native_value(vertex);
      }
      return best + offset(f, m);
    }
    const int j0 = std::max(0, static_cast<int>(std::floor((lo - grid_start) / grid_step)));
    const int j1 = std::min(grid_size - 1, static_cast<int>(std::ceil((hi - grid_start) / grid_step)));
    const double* values = &grid[(static_cast<size_t>(f) * groups + m) * grid_size];
    double best = R_NegInf;
    for (int j = j0; j <= j1; ++j) best = std::max(best, values[j]);
    return best;
  }
};

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
  for (int f = 0; f < groups; ++f)
    for (int m = 0; m < groups; ++m)
      for (int j = 0; j < draw.grid_size; ++j) {
        const size_t i = (static_cast<size_t>(f) * groups + m) * draw.grid_size + j;
        gaps[i] = draw.grid_start + j * draw.grid_step;
        female[i] = f;
        male[i] = m;
      }
  const NumericVector values = matching_utilities(gaps, model, -1, male, female);
  draw.grid.assign(values.begin(), values.end());
  return draw;
}

struct DrawWorkspace {
  std::vector<int> group, bin;
  std::vector<double> top, cumulative;
};

// Draw one partner for a woman of group `female_group` born `female_born`.
// Returns a person index, or -1 when no eligible man exists. `log_weights`
// holds one log group weight per group (zeros for female choice).
int draw_partner(const MalePool& pool, const DrawModel& draw, int female_group, double female_born,
                 double now, double cutoff, const std::vector<double>& log_weights,
                 const std::function<double(int)>& born_of,
                 const std::function<bool(int)>& excluded, DrawWorkspace& ws,
                 double& approximate) {
  if (pool.empty()) return -1;
  const int last = pool.last_bin(now, cutoff);
  const double w = pool.width();
  ws.group.clear(); ws.bin.clear(); ws.top.clear(); ws.cumulative.clear();
  double peak = R_NegInf;
  for (int g = 0; g < pool.groups(); ++g) {
    if (!draw.allowed(female_group, g)) continue;
    for (int k = 0; k <= last; ++k) {
      const std::vector<int>& members = pool.members(g, k);
      if (members.empty()) continue;
      const double oldest_birth = static_cast<double>(pool.first() + k) * w;
      const double top = draw.bound(female_group, g, female_born - oldest_birth - w,
                                    female_born - oldest_birth);
      const double value = log_weights[g] + top + std::log(static_cast<double>(members.size()));
      ws.group.push_back(g); ws.bin.push_back(k); ws.top.push_back(top);
      ws.cumulative.push_back(value);
      peak = std::max(peak, value);
    }
  }
  if (ws.group.empty() || peak == R_NegInf) return -1;
  double total = 0;
  for (double& value : ws.cumulative) { total += std::exp(value - peak); value = total; }
  for (int attempt = 0; attempt < 200; ++attempt) {
    const double u = R::runif(0, 1) * total;
    const size_t c = std::min(static_cast<size_t>(std::upper_bound(ws.cumulative.begin(),
      ws.cumulative.end(), u) - ws.cumulative.begin()), ws.group.size() - 1);
    const std::vector<int>& members = pool.members(ws.group[c], ws.bin[c]);
    const int man = members[std::min(static_cast<int>(members.size()) - 1,
      static_cast<int>(R::runif(0, 1) * members.size()))];
    const double utility = draw.utilities(female_group, ws.group[c],
      NumericVector::create(female_born - born_of(man)))[0];
    const double ratio = std::exp(utility - ws.top[c]);
    if (ratio > 1 + 1e-9) approximate += 1; // custom-formula bound was too low
    if (R::runif(0, 1) >= ratio) continue;
    if (now - born_of(man) < cutoff || excluded(man)) continue;
    return man;
  }
  // Many rejections (for example, most nearby men are kin): enumerate the
  // eligible pool and draw exactly from it.
  std::vector<int> candidates;
  std::vector<double> scores;
  for (int g = 0; g < pool.groups(); ++g) {
    if (!draw.allowed(female_group, g)) continue;
    std::vector<int> men;
    for (int k = 0; k <= last; ++k)
      for (int man : pool.members(g, k))
        if (now - born_of(man) >= cutoff && !excluded(man)) men.push_back(man);
    if (men.empty()) continue;
    NumericVector gaps(men.size());
    for (size_t j = 0; j < men.size(); ++j) gaps[j] = female_born - born_of(men[j]);
    const NumericVector utility = draw.utilities(female_group, g, gaps);
    for (size_t j = 0; j < men.size(); ++j) {
      candidates.push_back(men[j]);
      scores.push_back(log_weights[g] + utility[j]);
    }
  }
  if (candidates.empty()) return -1;
  const double best = *std::max_element(scores.begin(), scores.end());
  double sum = 0;
  for (double& score : scores) { sum += std::exp(score - best); score = sum; }
  const double u = R::runif(0, 1) * sum;
  const size_t j = std::min(static_cast<size_t>(std::upper_bound(scores.begin(), scores.end(), u) -
    scores.begin()), candidates.size() - 1);
  return candidates[j];
}

// Internal entry point: exact draws from a frozen pool, with the exact
// full-pool conditional-logit probabilities for comparison.
// [[Rcpp::export]]
List inspect_exact_draw_cpp(NumericVector male_born, IntegerVector male_group,
                            LogicalVector male_excluded, double female_born, int female_group,
                            NumericVector log_weights, List matching, int groups, double now,
                            double cutoff, int draws) {
  if (male_group.size() != male_born.size() || male_excluded.size() != male_born.size() ||
      log_weights.size() != groups || draws < 0)
    stop("Invalid exact-draw inputs.");
  double lo = std::min(female_born, now), hi = now;
  for (double born : male_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  MalePool pool(lo, hi, groups, kWeightBinWidth);
  for (R_xlen_t i = 0; i < male_born.size(); ++i) pool.add(i, male_group[i], male_born[i]);
  const DrawModel draw = make_draw_model(matching, groups, hi - lo);
  const std::vector<double> weights(log_weights.begin(), log_weights.end());
  auto born_of = [&](int man) { return male_born[man]; };
  auto excluded = [&](int man) { return static_cast<bool>(male_excluded[man]); };
  DrawWorkspace ws;
  double approximate = 0;
  IntegerVector counts(male_born.size());
  int failures = 0;
  for (int d = 0; d < draws; ++d) {
    if (d % 4096 == 0) checkUserInterrupt();
    const int man = draw_partner(pool, draw, female_group, female_born, now, cutoff, weights,
                                 born_of, excluded, ws, approximate);
    if (man < 0) ++failures; else ++counts[man];
  }
  NumericVector score(male_born.size(), R_NegInf);
  double best = R_NegInf;
  for (R_xlen_t i = 0; i < male_born.size(); ++i) {
    const int g = male_group[i];
    if (now - male_born[i] < cutoff || male_excluded[i] || !draw.allowed(female_group, g)) continue;
    score[i] = weights[g] + draw.utilities(female_group, g,
      NumericVector::create(female_born - male_born[i]))[0];
    best = std::max(best, score[i]);
  }
  NumericVector probability(male_born.size());
  double total = 0;
  for (R_xlen_t i = 0; i < male_born.size(); ++i)
    if (score[i] != R_NegInf) { probability[i] = std::exp(score[i] - best); total += probability[i]; }
  if (total > 0) for (double& p : probability) p /= total;
  return List::create(_["counts"] = counts, _["probability"] = probability,
                      _["failures"] = failures, _["approximate"] = approximate);
}

// Internal entry point: apply a sequence of male pool changes at increasing
// times, then compare incremental pool sums and weights with a full rebuild.
// [[Rcpp::export]]
List inspect_weighted_incremental_cpp(NumericVector female_born, IntegerVector female_groups,
                                      NumericVector male_born, IntegerVector male_groups,
                                      NumericVector op_times, IntegerVector op_men,
                                      NumericVector op_deltas, DataFrame rates, List matching,
                                      int group_count, double cutoff) {
  if (op_times.size() == 0 || op_men.size() != op_times.size() || op_deltas.size() != op_times.size())
    stop("Invalid incremental-calibration inputs.");
  double lo = op_times[0], hi = op_times[op_times.size() - 1];
  for (double born : female_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  for (double born : male_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  WeightCalibrator incremental(lo, hi, group_count, kWeightBinWidth, cutoff);
  const LagTable table = build_lag_table(matching, group_count, incremental.bins(), kWeightBinWidth);
  incremental.add_table(&table);
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    incremental.add_woman(female_groups[i], female_born[i], 1);
  }
  for (R_xlen_t i = 0; i < male_born.size(); ++i) incremental.add_man(male_groups[i], male_born[i], 1);
  std::vector<std::vector<Band>> bands(group_count);
  std::vector<const std::vector<Band>*> hazards(group_count);
  for (int g = 0; g < group_count; ++g) { bands[g] = read_bands(rates, "", "", g); hazards[g] = &bands[g]; }
  const std::vector<double> delta = as<std::vector<double>>(matching["delta"]);
  std::vector<double> warm;
  std::vector<double> prior(female_born.size(), 0);
  auto update_intensities = [&](double now) {
    for (R_xlen_t j = 0; j < female_born.size(); ++j) {
      const double value = current_hazard(hazards[female_groups[j]], female_born[j], now).first;
      incremental.replace_intensity(female_groups[j], female_born[j], {prior[j], prior[j]}, {value, value});
      prior[j] = value;
    }
  };
  for (R_xlen_t i = 0; i < op_times.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    update_intensities(op_times[i]);
    incremental.solve(op_times[i], false, 0, delta, warm);
    incremental.add_man(male_groups[op_men[i]], male_born[op_men[i]], op_deltas[i]);
  }
  const double end = op_times[op_times.size() - 1];
  WeightCalibrator rebuilt = incremental;
  rebuilt.refresh();
  const double difference = incremental.max_sum_difference(rebuilt, 0);
  std::vector<double> warm_a, warm_b;
  const WeightSolution a = incremental.solve(end, false, 0, delta, warm_a);
  const WeightSolution b = rebuilt.solve(end, false, 0, delta, warm_b);
  return List::create(_["max_sum_difference"] = difference,
                      _["incremental"] = wrap(a.log_weights), _["rebuilt"] = wrap(b.log_weights));
}

struct UnionRecord {
  double id, start, end;
  int wife, husband;
  bool wife_died = false;
  bool dissolved = false;
  bool dissolution_sampled = false;
  bool cohab = false, converted = false, conversion_sampled = false;
  double converted_from = NA_REAL, conversion_time = R_PosInf;
  double dissolution_time = std::numeric_limits<double>::infinity();
};

// [[Rcpp::export]]
List simulate_cpp(DataFrame persons, DataFrame mortality_rates,
                  DataFrame fertility_rates, double start, double years,
                  double birth_interval, DataFrame unions,
                  List matching, double min_male_age,
                  int kinship_generations, bool rate_initiation,
                  DataFrame union_rates, List father_matching,
                  bool cancel_birth_if_no_match, bool diagnostics,
                  NumericMatrix initial_ancestry, List inheritance, double ffab,
                  DataFrame dissolution_rates, DataFrame cohab_prob, DataFrame conversion_rates,
                  List cohab_matching, bool counts_cohab, bool check_cells = false,
                  Nullable<List> hetfert = R_NilValue, Nullable<DataFrame> transition_input = R_NilValue,
                  Nullable<List> model_schedule_input = R_NilValue) {
  const bool scheduled_models = model_schedule_input.isNotNull();
  List matching_bank, inheritance_bank;
  NumericVector model_times;
  IntegerVector general_route, father_route, cohab_route, inheritance_route;
  int model_epoch = 0, general_id = 0, father_id = 0, cohab_id = 0;
  if (scheduled_models) {
    List config(model_schedule_input.get());
    matching_bank = config["matching_models"];
    inheritance_bank = config["inheritance_models"];
    DataFrame epochs = as<DataFrame>(config["epochs"]);
    model_times = epochs["time_start"];
    general_route = epochs["matching"]; father_route = epochs["father_matching"];
    cohab_route = epochs["cohab_matching"]; inheritance_route = epochs["inheritance"];
    general_id = general_route[0]; father_id = father_route[0]; cohab_id = cohab_route[0];
  }
  // Keep all epochs in call-local tables; the population and its histories
  // stay in this single engine invocation across every rate change.
  const DataFrame all_mortality = mortality_rates, all_fertility = fertility_rates,
    all_search = union_rates, all_dissolution = dissolution_rates,
    all_conversion = conversion_rates, all_probability = cohab_prob;
  const DataFrame all_transitions = transition_input.isNotNull() ?
    DataFrame(transition_input.get()) : DataFrame::create();
  DataFrame active_transitions = calendar_slice(all_transitions, start);
  std::set<double> calendar_boundaries;
  for (const DataFrame& table : {all_mortality, all_fertility, all_search,
      all_dissolution, all_conversion, all_transitions, all_probability}) {
    if (!table.containsElementNamed("time_start")) continue;
    for (const char* column : {"time_start", "time_end"}) {
      NumericVector times = table[column];
      for (double time : times)
        if (start < time && time < start + years) calendar_boundaries.insert(time);
    }
  }
  if (scheduled_models) for (double time : model_times)
    if (start < time && time < start + years) calendar_boundaries.insert(time);
  mortality_rates = calendar_slice(all_mortality, start);
  fertility_rates = calendar_slice(all_fertility, start);
  union_rates = calendar_slice(all_search, start);
  dissolution_rates = calendar_slice(all_dissolution, start);
  conversion_rates = calendar_slice(all_conversion, start);
  cohab_prob = calendar_slice(all_probability, start);
  const List heterogeneity = hetfert.isNotNull() ? List(hetfert.get()) : List::create();
  const bool heterogeneous = heterogeneity.size() > 0;
  const double fertility_cv = heterogeneous ? as<double>(heterogeneity["cv"]) : 0;
  const double fertility_alpha = heterogeneous ? as<double>(heterogeneity["alpha"]) : 0;
  const double fertility_variance = fertility_cv <= 1 ? std::log1p(fertility_cv * fertility_cv) :
    2 * std::log(fertility_cv) + std::log1p((1 / fertility_cv) * (1 / fertility_cv));
  std::vector<double> fertility_weights;
  if (heterogeneous) {
    NumericVector supplied = persons["hetfert"];
    fertility_weights.assign(supplied.begin(), supplied.end());
  }
  const bool cohab_enabled = cohab_prob.nrows() > 0;
  const bool parity_specific = fertility_rates.containsElementNamed("parity_start");
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
  // Group indices cover starting identities and every declared possible one.
  int group_count = std::max(1, as<int>(inheritance["group_count"]));
  for (int value : person_groups) group_count = std::max(group_count, value + 1);
  // people retains every history, including deaths. pool holds living,
  // unpartnered men by group and birth bin for exact partner draws. index maps
  // public IDs to positions while reconstructing initial maternal histories.
  std::vector<Person> people;
  std::map<double, int> index;
  double next_id = 0, min_birth = start;
  for (int i = 0; i < ids.size(); ++i) {
    bool female = as<std::string>(sexes[i]) == "f";
    people.push_back({ids[i], births[i], deaths[i], moms[i], dads[i],
                      -std::numeric_limits<double>::infinity(), female, person_groups[i]});
    index[ids[i]] = i;
    next_id = std::max(next_id, ids[i]);
    min_birth = std::min(min_birth, births[i]);
  }
  std::vector<int> parity(parity_specific ? people.size() : 0, 0);
  MalePool pool(min_birth, start + years, group_count, kWeightBinWidth);
  for (int i = 0; i < static_cast<int>(people.size()); ++i)
    if (!people[i].female && NumericVector::is_na(people[i].died))
      pool.add(i, people[i].group, people[i].born);
  // Draw models for the union and (if separate) father matching models.
  // An empty population never draws; skip building custom-formula grids then.
  const DrawModel union_draw = people.empty() || scheduled_models ? DrawModel() :
    make_draw_model(matching, group_count, start + years - min_birth);
  const DrawModel* union_draw_current = &union_draw;
  std::map<int, DrawModel> scheduled_draws;
  auto get_scheduled_draw = [&](int id) -> const DrawModel* {
    auto found = scheduled_draws.find(id);
    if (found == scheduled_draws.end())
      found = scheduled_draws.emplace(id, people.empty() ? DrawModel() :
        make_draw_model(as<List>(matching_bank[id]), group_count, start + years - min_birth)).first;
    return &found->second;
  };
  DrawModel father_draw_storage;
  const DrawModel* father_draw = &union_draw;
  if (!scheduled_models && !people.empty() && SEXP(father_matching) != SEXP(matching)) {
    father_draw_storage = make_draw_model(father_matching, group_count, start + years - min_birth);
    father_draw = &father_draw_storage;
  }
  DrawModel cohab_draw_storage;
  const DrawModel* cohab_draw = &union_draw;
  if (!scheduled_models && cohab_enabled && !people.empty() && SEXP(cohab_matching) != SEXP(matching)) {
    cohab_draw_storage = make_draw_model(cohab_matching, group_count, start + years - min_birth);
    cohab_draw = &cohab_draw_storage;
  }
  if (scheduled_models) {
    union_draw_current = get_scheduled_draw(general_id);
    father_draw = get_scheduled_draw(father_id);
    cohab_draw = get_scheduled_draw(cohab_id);
  }
  DrawWorkspace draw_workspace;
  double draw_approximate = 0;
  std::vector<std::array<int, 2>> parent_links(people.size());
  for (int person = 0; person < static_cast<int>(people.size()); ++person) {
    parent_links[person] = {
      NumericVector::is_na(people[person].mother) ? -1 : index.at(people[person].mother),
      NumericVector::is_na(people[person].father) ? -1 : index.at(people[person].father)
    };
  }
  hammel::Kinship kinship(std::move(parent_links), kinship_generations);
  std::vector<UnionRecord> union_history;
  // Log group weights applied at each union-forming search (empty = unit weights).
  std::vector<std::vector<double>> diagnostic_log_weights;
  std::vector<std::vector<int>> person_unions(people.size());
  std::vector<int> current_partner(people.size(), -1);
  const bool previous_specific = union_rates.containsElementNamed("previous_union");
  std::vector<int> previous_union(people.size(), 0), legal_previous(people.size(), 0);
  std::vector<double> previous_union_end(people.size(), NA_REAL);
  if (previous_specific || cohab_enabled) {
    IntegerVector initial_previous = persons[".hammel_previous_union"];
    previous_union = as<std::vector<int>>(initial_previous);
    previous_union_end = as<std::vector<double>>(persons[".hammel_previous_union_end"]);
    legal_previous = as<std::vector<int>>(persons[".hammel_legal_previous"]);
  }
  // weighted_choice state (D066). The calibrator exists only for grouped
  // populations; until it is created below, pool hooks record nothing.
  std::unique_ptr<WeightCalibrator> calibrator;
  std::vector<char> counted_woman;
  std::unique_ptr<FemaleIntensityTracker> intensity_tracker;
  // Cell-accounted fertility under birth spacing (D101).
  std::unique_ptr<FertilityCells> cells;
  std::vector<const std::vector<Band>*> search_hazards;
  std::vector<bool> search_duration;
  std::vector<std::vector<Band>> probability_schedules(group_count * 3);
  auto load_probabilities = [&]() {
    if (cohab_enabled) {
      const bool previous_prob = cohab_prob.containsElementNamed("previous_union");
      const std::array<std::string, 3> labels = {"none", "separated", "widowed"};
      for (int g = 0; g < group_count; ++g) for (int p = 0; p < 3; ++p)
        probability_schedules[g * 3 + p] = read_bands(cohab_prob, "", "", g, previous_prob ? labels[p] : "");
    }
  };
  load_probabilities();
  // Pool updates are constant-time (swap removal within a cell).
  auto remove_male = [&](int man) {
    if (!pool.contains(man)) return; // Partnered men are already absent when they die.
    pool.remove(man);
    if (calibrator) calibrator->add_man(people[man].group, people[man].born, -1);
  };
  auto add_male = [&](int man) {
    if (people[man].female || !NumericVector::is_na(people[man].died) ||
        current_partner[man] >= 0 || pool.contains(man)) return;
    pool.add(man, people[man].group, people[man].born);
    if (calibrator) calibrator->add_man(people[man].group, people[man].born, 1);
  };
  double next_union_id = 0;
  NumericVector uids = unions["union_id"], wives = unions["wife_id"];
  NumericVector husbands = unions["husband_id"], starts = unions["start_time"];
  NumericVector ends = unions["end_time"];
  for (int u = 0; u < uids.size(); ++u) {
    const int woman = index.at(wives[u]), man = index.at(husbands[u]);
    union_history.push_back({uids[u], starts[u], ends[u], woman, man});
    if (unions.containsElementNamed("union_type")) {
      CharacterVector types = unions["union_type"];
      union_history.back().cohab = as<std::string>(types[u]) == "cohab";
    }
    if (unions.containsElementNamed("previous_cohab_id")) {
      NumericVector links = unions["previous_cohab_id"];
      union_history.back().converted_from = links[u];
    }
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
      if (found != index.end()) {
        people[found->second].last_birth = std::max(people[found->second].last_birth, child.born);
        if (parity_specific && people[found->second].female) ++parity[found->second];
      }
    }
  }
  // weighted_choice counts living, unpartnered women; fertility-driven intensity
  // also requires the birth-spacing interval to have ended (queued re-entry).
  auto add_counted_woman = [&](int woman, double now) {
    if (!calibrator || !people[woman].female || !NumericVector::is_na(people[woman].died) ||
        current_partner[woman] >= 0 || counted_woman[woman]) return;
    counted_woman[woman] = 1;
    calibrator->add_woman(people[woman].group, people[woman].born, 1);
    const int g = people[woman].group;
    const int index = previous_specific ? g * 3 + previous_union[woman] : g;
    const auto* search = !rate_initiation ? nullptr : search_hazards[index];
    const double origin = rate_initiation && search_duration[index] ?
      previous_union_end[woman] : people[woman].born;
    if (!R_finite(origin)) stop("Duration search risk requires a recorded previous union ending.");
    // The fertility channel comes from cells (fill_calibration), not women.
    intensity_tracker->activate(woman, g, people[woman].born, search, nullptr,
      people[woman].last_birth + birth_interval, now, origin,
      cohab_enabled ? &probability_schedules[g * 3 + previous_union[woman]] : nullptr);
  };
  auto remove_counted_woman = [&](int woman) {
    if (!calibrator || !counted_woman[woman]) return;
    counted_woman[woman] = 0;
    intensity_tracker->remove(woman);
    calibrator->add_woman(people[woman].group, people[woman].born, -1);
  };
  auto counted_birth = [&](int woman, double now) {
    if (calibrator && counted_woman[woman]) {
      if (parity_specific) { remove_counted_woman(woman); add_counted_woman(woman, now); }
      else intensity_tracker->set_ready(woman, now + birth_interval, now);
    }
  };
  LagTable union_table, father_table_storage, cohab_table;
  int union_table_id = 0, father_table_id = 0, cohab_table_id = 0;
  // A NULL father model shares the general model and calibration table. When a
  // model prevents all exogamy, every woman sees only her own group's men, so
  // a group weight cannot change her choice: skip calibration.
  auto applies_weights = [](const List& model) {
    return !as<bool>(model["prevent_exogamy_all"]);
  };
  bool union_weighted = applies_weights(matching);
  bool father_weighted = applies_weights(father_matching);
  bool cohab_weighted = cohab_enabled && applies_weights(cohab_matching);
  bool any_weighted = union_weighted || father_weighted || cohab_weighted;
  if (scheduled_models) for (auto model : matching_bank)
    any_weighted = any_weighted || applies_weights(as<List>(model));
  std::map<int, LagTable> scheduled_lags;
  auto get_scheduled_lag = [&](int id) -> const LagTable* {
    auto found = scheduled_lags.find(id);
    if (found == scheduled_lags.end())
      found = scheduled_lags.emplace(id, build_lag_table(as<List>(matching_bank[id]),
        group_count, calibrator->bins(), kWeightBinWidth)).first;
    return &found->second;
  };
  std::vector<double> union_warm, father_warm, cohab_warm;
  double weighted_infeasible = 0;
  // Partner and father searches share eligibility and the exact draw; only
  // form_union changes partner state or removes the selected man.
  auto find_candidate = [&](int woman, double now, const List& choice_model,
                            const DrawModel& draw, bool father,
                            std::vector<double>* applied, bool cohab_role = false) -> int {
    if (pool.empty()) return -1;
    std::vector<double> log_weights(group_count, 0.0);
    if (calibrator && (cohab_role ? cohab_weighted : (father ? father_weighted : union_weighted))) {
      // Searches and on-demand unions use the union-formation calibration;
      // rate-mode nonmarital fathers have their own, driven by fertility.
      intensity_tracker->advance(now);
      const bool fertility_driven = !rate_initiation || father;
      if (cells && !cohab_role && fertility_driven) {
        const std::vector<double> filled = cells->fill_calibration(*calibrator);
        if (check_cells) {
          // Explicit scan of counted women's individual hazards.
          std::vector<double> scan(filled.size(), 0.0);
          for (int w = 0; w < static_cast<int>(counted_woman.size()); ++w)
            if (counted_woman[w]) scan[calibrator->bin_index(people[w].group, people[w].born)] += cells->hazard(w);
          for (size_t k = 0; k < scan.size(); ++k)
            if (std::abs(scan[k] - filled[k]) > 1e-9 * std::max(1.0, std::abs(scan[k])))
              stop("Internal error: cell calibration totals differ from an explicit scan.");
        }
      }
      const WeightSolution weights = calibrator->solve(now, cohab_role ? 2 : (fertility_driven ? 1 : 0),
        cohab_role ? cohab_table_id : (father ? father_table_id : union_table_id),
        as<std::vector<double>>(choice_model["delta"]), cohab_role ? cohab_warm : (father ? father_warm : union_warm));
      if (weights.infeasible) ++weighted_infeasible;
      log_weights = weights.log_weights;
    }
    if (applied) *applied = log_weights;
    return draw_partner(pool, draw, people[woman].group, people[woman].born, now, min_male_age,
      log_weights, [&](int person) { return people[person].born; },
      [&](int man) { return kinship.generations() > 0 && kinship.known_related(woman, man); },
      draw_workspace, draw_approximate);
  };
  auto form_union = [&](int woman, double now) {
    int partner = current_partner[woman];
    if (partner >= 0) return partner;
    std::vector<double> applied;
    const double p = cohab_enabled ? current_hazard(
      &probability_schedules[people[woman].group * 3 + previous_union[woman]], people[woman].born, now).first : 0;
    const bool is_cohab = cohab_enabled && R::runif(0, 1) < p;
    partner = find_candidate(woman, now, is_cohab ? cohab_matching : matching,
      is_cohab ? *cohab_draw : *union_draw_current, false, diagnostics ? &applied : nullptr, is_cohab);
    if (partner >= 0) {
      if (next_union_id >= 9007199254740991.0 ||
          union_history.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Union history exceeds supported ID or row limits.");
      const int u = static_cast<int>(union_history.size());
      union_history.push_back({++next_union_id, now, NA_REAL, woman, partner});
      union_history.back().cohab = is_cohab;
      if (diagnostics) {
        diagnostic_log_weights.push_back(applied);
      }
      person_unions[woman].push_back(u);
      person_unions[partner].push_back(u);
      current_partner[woman] = partner;
      current_partner[partner] = woman;
      remove_male(partner);
      remove_counted_woman(woman);
    }
    return partner;
  };
  struct GroupSchedules {
    std::vector<Band> male, female, search, conversion;
    bool duration_conversion = false;
    std::array<std::vector<Band>, 3> previous_search;
    std::array<bool, 3> duration_search = {false, false, false};
  };
  int inheritance_method = inheritance["method"];
  double passing = inheritance["passing"];
  IntegerVector inheritance_ranks = inheritance["ranks"];
  IntegerVector ancestry_groups = inheritance["ancestry_groups"];
  Function inherit_custom = inheritance["custom"];
  std::vector<GroupSchedules> schedules(group_count);
  const bool status_fertility = fertility_rates.containsElementNamed("union_type");
  bool typed_fertility = false;
  auto load_schedules = [&]() {
    for (int group = 0; group < group_count; ++group) {
      auto& rates = schedules[group];
      rates.male = read_bands(mortality_rates, "m", "", group);
      rates.female = read_bands(mortality_rates, "f", "", group);
      if (conversion_rates.nrows()) {
        NumericVector ages = conversion_rates["age_start"];
        rates.duration_conversion = NumericVector::is_na(ages[0]);
        rates.conversion = read_bands(conversion_rates, "", "", group, "", rates.duration_conversion);
      }
      if (rate_initiation) {
        if (previous_specific) {
          const std::array<std::string, 3> labels = {"none", "separated", "widowed"};
          for (int p = 0; p < 3; ++p) {
            rates.duration_search[p] = search_uses_duration(union_rates, labels[p]);
            rates.previous_search[p] = read_bands(union_rates, "", "", group, labels[p],
              rates.duration_search[p]);
          }
        } else rates.search = read_bands(union_rates, "", "", group);
      }
    }
  };
  load_schedules();
  // weighted_choice needs schedules and group indices, so it starts here.
  // One group means unit weights, so calibration is unnecessary.
  if (any_weighted && group_count > 1) {
    double min_birth = start;
    for (const Person& person : people) min_birth = std::min(min_birth, person.born);
    calibrator.reset(new WeightCalibrator(min_birth, start + years, group_count,
                                          kWeightBinWidth, min_male_age));
    for (int group = 0; group < group_count; ++group) {
      if (previous_specific) {
        for (int p = 0; p < 3; ++p) {
          search_hazards.push_back(&schedules[group].previous_search[p]);
          search_duration.push_back(schedules[group].duration_search[p]);
        }
      } else {
        search_hazards.push_back(&schedules[group].search);
        search_duration.push_back(false);
      }
    }
    if (scheduled_models) {
      union_table_id = calibrator->add_table(get_scheduled_lag(general_id));
      father_table_id = calibrator->add_table(get_scheduled_lag(father_id));
      cohab_table_id = calibrator->add_table(get_scheduled_lag(cohab_id));
    } else {
      if (union_weighted) {
        union_table = build_lag_table(matching, group_count, calibrator->bins(), kWeightBinWidth);
        union_table_id = father_table_id = calibrator->add_table(&union_table);
      }
      if (father_weighted && SEXP(father_matching) != SEXP(matching)) {
        father_table_storage = build_lag_table(father_matching, group_count, calibrator->bins(),
                                               kWeightBinWidth);
        father_table_id = calibrator->add_table(&father_table_storage);
      }
      if (cohab_weighted) {
        if (SEXP(cohab_matching) == SEXP(matching)) cohab_table_id = union_table_id;
        else {
          cohab_table = build_lag_table(cohab_matching, group_count, calibrator->bins(), kWeightBinWidth);
          cohab_table_id = calibrator->add_table(&cohab_table);
        }
      }
    }
    counted_woman.assign(people.size(), 0);
    intensity_tracker.reset(new FemaleIntensityTracker(*calibrator));
    for (int man = 0; man < static_cast<int>(people.size()); ++man)
      if (pool.contains(man)) calibrator->add_man(people[man].group, people[man].born, 1);
    for (int woman = 0; woman < static_cast<int>(people.size()); ++woman)
      add_counted_woman(woman, start);
  }
  auto update_models = [&](double now) {
    if (!scheduled_models) return;
    while (model_epoch + 1 < model_times.size() && model_times[model_epoch + 1] <= now)
      ++model_epoch;
    const bool matching_changed = general_id != general_route[model_epoch] ||
      father_id != father_route[model_epoch] || cohab_id != cohab_route[model_epoch];
    if (matching_changed) {
      general_id = general_route[model_epoch]; father_id = father_route[model_epoch];
      cohab_id = cohab_route[model_epoch];
      matching = as<List>(matching_bank[general_id]);
      father_matching = as<List>(matching_bank[father_id]);
      cohab_matching = as<List>(matching_bank[cohab_id]);
      union_draw_current = get_scheduled_draw(general_id);
      father_draw = get_scheduled_draw(father_id); cohab_draw = get_scheduled_draw(cohab_id);
      union_weighted = applies_weights(matching); father_weighted = applies_weights(father_matching);
      cohab_weighted = cohab_enabled && applies_weights(cohab_matching);
      union_warm.clear(); father_warm.clear(); cohab_warm.clear();
      if (calibrator) {
        calibrator->set_table(union_table_id, get_scheduled_lag(general_id));
        calibrator->set_table(father_table_id, get_scheduled_lag(father_id));
        calibrator->set_table(cohab_table_id, get_scheduled_lag(cohab_id));
        calibrator->refresh();
      }
    }
    inheritance = as<List>(inheritance_bank[inheritance_route[model_epoch]]);
    inheritance_method = as<int>(inheritance["method"]);
    passing = as<double>(inheritance["passing"]);
    inheritance_ranks = as<IntegerVector>(inheritance["ranks"]);
    ancestry_groups = as<IntegerVector>(inheritance["ancestry_groups"]);
    inherit_custom = as<Function>(inheritance["custom"]);
  };
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
  const bool has_dissolution = dissolution_rates.nrows() > 0;
  std::vector<std::vector<Band>> dissolution_schedules, cohab_dissolution_schedules;
  auto load_dissolution = [&]() {
    dissolution_schedules.clear(); cohab_dissolution_schedules.clear();
    if (has_dissolution) {
      const int count = dissolution_rates.containsElementNamed(".hammel_group") ?
        group_count * group_count : 1;
      const bool typed = dissolution_rates.containsElementNamed("union_type");
      for (int pair = 0; pair < count; ++pair) {
        dissolution_schedules.push_back(read_bands(dissolution_rates, "", typed ? "marriage" : "", pair));
        cohab_dissolution_schedules.push_back(read_bands(dissolution_rates, "", typed ? "cohab" : "", pair));
      }
    }
  };
  load_dissolution();
  struct TransitionSchedule { int dest, sex; std::string type; std::vector<Band> bands; };
  std::vector<std::vector<TransitionSchedule>> transition_schedules(group_count);
  bool transition_duration = false;
  const bool has_transition = transition_input.isNotNull();
  auto load_transitions = [&]() {
    for (auto& origin : transition_schedules) origin.clear();
    if (has_transition) {
      DataFrame input = active_transitions;
      transition_duration = input.containsElementNamed("dur_start");
      IntegerVector origins = input[".origin"], destinations = input[".dest"];
      NumericVector lo = input[transition_duration ? "dur_start" : "age_start"];
      NumericVector hi = input[transition_duration ? "dur_end" : "age_end"], values = input["rate"];
      const bool sexed = input.containsElementNamed("sex"), typed = input.containsElementNamed("union_type");
      CharacterVector sexes = sexed ? as<CharacterVector>(input["sex"]) : CharacterVector();
      CharacterVector types = typed ? as<CharacterVector>(input["union_type"]) : CharacterVector();
      for (int row = 0; row < input.nrows(); ++row) {
        auto& schedules = transition_schedules[origins[row]];
        const int sex = !sexed ? -1 : (as<std::string>(sexes[row]) == "f" ? 1 : 0);
        const std::string type = typed ? as<std::string>(types[row]) : "";
        auto found = std::find_if(schedules.begin(), schedules.end(), [&](const TransitionSchedule& x) {
          return x.dest == destinations[row] && x.sex == sex && x.type == type;
        });
        if (found == schedules.end()) {
          schedules.push_back({destinations[row], sex, type, {}});
          found = std::prev(schedules.end());
        }
        found->bands.push_back({lo[row], hi[row], values[row]});
      }
      for (auto& origin : transition_schedules) for (auto& schedule : origin)
        std::sort(schedule.bands.begin(), schedule.bands.end(), [](const Band& a, const Band& b) { return a.start < b.start; });
    }
  };
  load_transitions();
  std::vector<double> group_entry(people.size());
  if (has_transition) group_entry = as<std::vector<double>>(persons[".group_entry"]);
  std::vector<int> transition_destination(people.size(), -1);
  std::vector<double> tx_person, tx_time;
  std::vector<int> tx_origin, tx_dest;
  auto schedule = [&](int i, double now) {
    const auto version = ++event_versions[i];
    const Person& p = people[i];
    if (p.female) cells->refresh(i, now);
    if (!NumericVector::is_na(p.died)) return;
    const auto& rates = schedules[p.group];
    const auto& death_rates = p.female ? rates.female : rates.male;
    double death_time = candidate(waiting_time(death_rates, now - p.born, threshold()), now);
    // Births come from fertility cells (D101), not per-person candidates.
    const double birth_time = std::numeric_limits<double>::infinity();
    double search_time = std::numeric_limits<double>::infinity();
    if (rate_initiation && p.female && current_partner[i] < 0) {
      const bool duration = previous_specific && rates.duration_search[previous_union[i]];
      const double origin = duration ? previous_union_end[i] : p.born;
      if (!R_finite(origin)) stop("Duration search risk requires a recorded previous union ending.");
      search_time = candidate(waiting_time(previous_specific ? rates.previous_search[previous_union[i]] : rates.search,
        now - origin, threshold()), now);
    }
    double dissolution_time = std::numeric_limits<double>::infinity();
    if (has_dissolution && p.female && current_partner[i] >= 0) {
      UnionRecord& relationship = union_history[person_unions[i].back()];
      if (!relationship.dissolution_sampled) {
        const int pair = dissolution_schedules.size() == 1 ? 0 :
          p.group + group_count * people[current_partner[i]].group;
        relationship.dissolution_time = candidate(waiting_time(relationship.cohab ? cohab_dissolution_schedules[pair] : dissolution_schedules[pair],
          now - relationship.start, threshold()), now);
        relationship.dissolution_sampled = true;
      }
      // Retain the union clock through births, including a dissolution due now
      // after a tied partnered birth. It is not a new zero-wait random draw.
      dissolution_time = relationship.dissolution_time;
    }
    double conversion_time = R_PosInf;
    if (cohab_enabled && !rates.conversion.empty() && p.female && current_partner[i] >= 0) {
      auto& relationship = union_history[person_unions[i].back()];
      if (relationship.cohab) {
        if (!relationship.conversion_sampled) {
          relationship.conversion_time = candidate(waiting_time(rates.conversion,
            now - (rates.duration_conversion ? relationship.start : p.born), threshold()), now);
          relationship.conversion_sampled = true;
        }
        conversion_time = relationship.conversion_time;
      }
    }
    double transition_time = R_PosInf;
    transition_destination[i] = -1;
    if (has_transition) {
      const int partner = current_partner[i];
      const std::string type = partner < 0 ? "unpartnered" :
        (union_history[person_unions[i].back()].cohab ? "cohab" : "marriage");
      for (const auto& outgoing : transition_schedules[p.group]) {
        if (outgoing.sex >= 0 && outgoing.sex != (p.female ? 1 : 0)) continue;
        if (!outgoing.type.empty() && outgoing.type != type &&
            !(outgoing.type == "partnered" && partner >= 0)) continue;
        const double due = candidate(waiting_time(outgoing.bands,
          now - (transition_duration ? group_entry[i] : p.born), threshold()), now);
        if (due < transition_time) {
          transition_time = due;
          transition_destination[i] = outgoing.dest;
        }
      }
    }
    const double other_due = std::min(std::min(death_time, birth_time),
      std::min(std::min(dissolution_time, conversion_time), search_time));
    const double due = std::min(other_due, transition_time);
    const auto type = death_time <= due ? EventType::death :
      (birth_time <= due ? EventType::birth :
       (dissolution_time <= due ? EventType::dissolution :
        (conversion_time <= due ? EventType::conversion :
         (search_time <= due ? EventType::search : EventType::transition))));
    if (due < end) queue.emplace(due, Event{i, type, version});
  };
  // Fertility type labels: overall, unpartnered/partnered, or three types.
  std::vector<std::string> fertility_types;
  auto load_cells = [&](double now) {
    typed_fertility = false;
    if (status_fertility) {
      CharacterVector types = fertility_rates["union_type"];
      for (auto type : types) if (as<std::string>(type) == "cohab") typed_fertility = true;
    }
    fertility_types = !status_fertility ? std::vector<std::string>{""} :
      (typed_fertility ? std::vector<std::string>{"unpartnered", "marriage", "cohab"} :
                         std::vector<std::string>{"unpartnered", "partnered"});
    auto state = [&](int w) {
      FertilityCells::State st;
      const Person& p = people[w];
      if (!p.female || !NumericVector::is_na(p.died)) return st;
      int type = 0;
      if (status_fertility && current_partner[w] >= 0)
        type = typed_fertility && union_history[person_unions[w].back()].cohab ? 2 : 1;
      st.stratum = cells->stratum(p.group, type, parity_specific ? parity[w] : 0);
      if (st.stratum < 0) stop("No fertility schedule covers a woman's current parity.");
      st.weight = heterogeneous ? fertility_weights[w] : 1;
      st.born = p.born;
      st.ready = p.last_birth + birth_interval;
      st.unpartnered = current_partner[w] < 0;
      return st;
    };
    double oldest_birth = end;
    for (const Person& person : people) oldest_birth = std::min(oldest_birth, person.born);
    cells.reset(new FertilityCells(fertility_strata(fertility_rates, group_count, fertility_types),
      group_count, static_cast<int>(fertility_types.size()), parity_specific, end - oldest_birth,
      birth_interval, heterogeneous, state, threshold, calibrator.get(), now, end));
  };
  load_cells(start);
  // Initialize one next event for each living person. Known age changes are
  // already incorporated, so birthdays alone do not trigger resampling.
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    schedule(i, start);
  }
  double birth_count = 0, death_count = 0, dissolution_count = 0, conversion_count = 0;
  unsigned long long processed = 0;
  std::vector<DataFrame> unrealized_epochs;
  bool cells_consistent = true;
  // Process events in chronological order until no scheduled events remain.
  while (true) {
    // Calendar changes win exact ties with every cell/person boundary or event.
    // Conditional Exp(1) resampling is valid after survival to the boundary.
    // Candidates under the superseded rates never execute beyond it.
    if (!calendar_boundaries.empty()) {
      const double now = *calendar_boundaries.begin();
      if (now <= cells->next_boundary() && now <= cells->next_clock() &&
          (queue.empty() || now <= queue.begin()->first)) {
        calendar_boundaries.erase(calendar_boundaries.begin());
        checkUserInterrupt();
        update_models(now);
        for (int w = 0; w < static_cast<int>(counted_woman.size()); ++w)
          remove_counted_woman(w);
        cohab_prob = calendar_slice(all_probability, now);
        load_probabilities();
        const bool rate_change = calendar_changes_at(all_mortality, now) ||
          calendar_changes_at(all_fertility, now) || calendar_changes_at(all_search, now) ||
          calendar_changes_at(all_dissolution, now) || calendar_changes_at(all_conversion, now) ||
          calendar_changes_at(all_transitions, now);
        if (!rate_change) {
          // Model/probability changes retain all demographic event clocks.
          // Role preferences and intensities use the newly active inputs.
          for (int w = 0; w < static_cast<int>(people.size()); ++w) {
            if (w % 1024 == 0) checkUserInterrupt();
            add_counted_woman(w, now);
          }
          continue;
        }
        const bool fertility_change = calendar_changes_at(all_fertility, now);
        if (fertility_change) {
          cells->finish(now);
          cells_consistent = cells_consistent && cells->consistent();
          unrealized_epochs.push_back(cells->unrealized(fertility_types, now));
        }
        mortality_rates = calendar_slice(all_mortality, now);
        fertility_rates = calendar_slice(all_fertility, now);
        union_rates = calendar_slice(all_search, now);
        dissolution_rates = calendar_slice(all_dissolution, now);
        conversion_rates = calendar_slice(all_conversion, now);
        active_transitions = calendar_slice(all_transitions, now);
        load_schedules(); load_dissolution(); load_transitions();
        if (fertility_change) load_cells(now);
        for (auto& relationship : union_history) {
          if (!NumericVector::is_na(relationship.end)) continue;
          relationship.dissolution_sampled = false;
          relationship.conversion_sampled = false;
        }
        // All table objects referenced by the tracker keep their addresses.
        // Clock selection is fixed across calendar periods by R validation.
        for (int w = 0; w < static_cast<int>(people.size()); ++w)
          add_counted_woman(w, now);
        queue.clear();
        for (int i = 0; i < static_cast<int>(people.size()); ++i) {
          if (i % 1024 == 0) checkUserInterrupt();
          schedule(i, now);
        }
        continue;
      }
    }
    // Cell clocks depend on current membership, so bookkeeping boundaries are
    // processed in time order before any event at or after them.
    if (cells) {
      const double boundary = cells->next_boundary();
      if (boundary < end && boundary <= cells->next_clock() &&
          (queue.empty() || boundary <= queue.begin()->first)) {
        cells->advance_boundary();
        continue;
      }
    }
    // Exact ties follow the event priority death, birth, dissolution,
    // conversion, search: person events at the same time as a cell clock go
    // first only while a current death is pending at that time.
    const double clock = cells ? cells->next_clock() : R_PosInf;
    if (queue.empty() && !(clock < end)) break;
    if (++processed % 1024 == 0) checkUserInterrupt();
    double now;
    EventType type;
    int i;
    auto death_at = [&](double time) {
      const auto range = queue.equal_range(time);
      for (auto it = range.first; it != range.second; ++it)
        if (it->second.type == EventType::death && it->second.version == event_versions[it->second.person])
          return true;
      return false;
    };
    if (clock < end && (queue.empty() || clock < queue.begin()->first ||
        (clock == queue.begin()->first && !death_at(clock)))) {
      now = clock;
      i = cells->fire_next(now);
      type = EventType::birth;
    } else {
      auto next = queue.begin();
      now = next->first;
      const Event event = next->second;
      queue.erase(next);
      type = event.type;
      i = event.person;
      if (event.version != event_versions[i]) continue;
    }
    if (intensity_tracker) intensity_tracker->advance(now);
    if (check_cells && cells && !cells->consistent())
      stop("Internal error: fertility cell counts differ from an explicit scan.");
    if (type == EventType::transition) {
      const int destination = transition_destination[i];
      if (destination < 0 || destination == people[i].group) stop("Invalid transition destination.");
      const int partner = current_partner[i];
      if (people[i].female) remove_counted_woman(i); else remove_male(i);
      tx_person.push_back(people[i].id); tx_time.push_back(now);
      tx_origin.push_back(people[i].group); tx_dest.push_back(destination);
      people[i].group = destination;
      group_entry[i] = now;
      if (people[i].female) add_counted_woman(i, now); else add_male(i);
      if (partner >= 0) {
        auto& relationship = union_history[person_unions[i].back()];
        relationship.dissolution_sampled = false;
        if (people[i].female) relationship.conversion_sampled = false;
        // The woman carries the pair-dependent union clock even when he moves.
        if (!people[i].female) schedule(partner, now);
      }
      schedule(i, now);
    } else if (type == EventType::death) {
      // Retain the person's history, record death, and stop their scheduling.
      people[i].died = now;
      for (int u : person_unions[i]) {
        UnionRecord& relationship = union_history[u];
        if (!NumericVector::is_na(relationship.end)) continue;
        relationship.end = now;
        relationship.wife_died = i == relationship.wife;
        current_partner[relationship.wife] = -1;
        current_partner[relationship.husband] = -1;
        if (!relationship.cohab) legal_previous[relationship.wife] = legal_previous[relationship.husband] = 2;
        previous_union[relationship.wife] = relationship.cohab && !counts_cohab ? legal_previous[relationship.wife] : 2;
        previous_union[relationship.husband] = relationship.cohab && !counts_cohab ? legal_previous[relationship.husband] : 2;
        previous_union_end[relationship.wife] = previous_union_end[relationship.husband] = now;
        add_male(relationship.husband);
        add_counted_woman(relationship.wife, now);
        if (has_transition && i == relationship.wife) schedule(relationship.husband, now);
        if (cells) cells->refresh(relationship.wife, now);
        // Her old birth candidate is invalid once her husband has died.
        if ((rate_initiation || has_dissolution || has_transition) && i == relationship.husband)
          schedule(relationship.wife, now);
      }
      if (!people[i].female) remove_male(i);
      else remove_counted_woman(i);
      if (cells && people[i].female) cells->refresh(i, now);
      ++death_count;
    } else if (type == EventType::dissolution) {
      const int partner = current_partner[i];
      if (partner < 0) stop("Dissolution requires an ongoing union.");
      UnionRecord& relationship = union_history[person_unions[i].back()];
      relationship.end = now;
      relationship.dissolved = true;
      current_partner[i] = current_partner[partner] = -1;
      if (!relationship.cohab) legal_previous[i] = legal_previous[partner] = 1;
      previous_union[i] = relationship.cohab && !counts_cohab ? legal_previous[i] : 1;
      previous_union[partner] = relationship.cohab && !counts_cohab ? legal_previous[partner] : 1;
      previous_union_end[i] = previous_union_end[partner] = now;
      add_male(partner);
      add_counted_woman(i, now);
      ++dissolution_count;
      // The woman's cell/search state changes. With transition risks enabled,
      // both partners also need fresh union-type-dependent candidates.
      schedule(i, now);
      if (has_transition) schedule(partner, now);
    } else if (type == EventType::conversion) {
      const int partner = current_partner[i];
      const int old = person_unions[i].back();
      if (partner < 0 || !union_history[old].cohab) stop("Conversion requires an ongoing cohabitation.");
      union_history[old].end = now;
      union_history[old].converted = true;
      const double prior_id = union_history[old].id;
      if (next_union_id >= 9007199254740991.0 || union_history.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Union history exceeds supported ID or row limits.");
      const int next = union_history.size();
      union_history.push_back({++next_union_id, now, NA_REAL, i, partner});
      union_history.back().converted_from = prior_id;
      person_unions[i].push_back(next);
      person_unions[partner].push_back(next);
      if (diagnostics) diagnostic_log_weights.emplace_back(group_count, NA_REAL);
      ++conversion_count;
      schedule(i, now);
      if (has_transition) schedule(partner, now);
    } else if (type == EventType::search) {
      const int partner = form_union(i, now);
      schedule(i, now);
      if (has_transition && partner >= 0) schedule(partner, now);
    } else {
      double father = NA_REAL;
      int partner = current_partner[i];
      if (!rate_initiation && partner < 0) {
        partner = form_union(i, now);
        if (has_transition && partner >= 0) schedule(partner, now);
      }
      // A nonmarital father stays available: parentage alone is not a union.
      const int father_position = partner >= 0 ? partner :
        (rate_initiation ? find_candidate(i, now, father_matching, *father_draw, true, nullptr) : -1);
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
      if (parity_specific) ++parity[i];
      counted_birth(i, now);
      // Default ffab is 0.4886 - its the signature of Berkeley training!
      bool female = R::runif(0, 1) < ffab;
      int child = static_cast<int>(people.size());
      people.push_back({++next_id, now, NA_REAL, mother, father,
                        -std::numeric_limits<double>::infinity(), female, people[i].group});
      if (parity_specific) parity.push_back(0);
      if (heterogeneous) fertility_weights.push_back(female ?
        inherit_fertility_weight(fertility_weights[i], fertility_variance, fertility_alpha) : NA_REAL);
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
      // hazard or matching lookup. Transitions separately update existing identities.
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
      if (!female) {
        pool.add(child, child_group, now);
        if (calibrator) calibrator->add_man(child_group, now, 1);
      }
      person_unions.emplace_back();
      current_partner.push_back(-1);
      previous_union.push_back(0);
      legal_previous.push_back(0);
      previous_union_end.push_back(NA_REAL);
      if (calibrator) {
        counted_woman.push_back(0);
        add_counted_woman(child, now);
      }
      event_versions.push_back(0);
      group_entry.push_back(now);
      transition_destination.push_back(-1);
      ++birth_count;
      // Both newborn and mother need fresh candidates. The mother now has a
      // new last_birth; spacing delays her next birth, but not a possible death.
      schedule(child, now);
      schedule(i, now);
    }
  }
  DataFrame unrealized;
  if (cells) {
    cells->finish(end);
    cells_consistent = cells_consistent && cells->consistent();
    unrealized = cells->unrealized(fertility_types);
    Function bind = Environment::base_env()["rbind"];
    for (auto it = unrealized_epochs.rbegin(); it != unrealized_epochs.rend(); ++it)
      unrealized = as<DataFrame>(bind(*it, unrealized));
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
  CharacterVector out_reason(union_history.size()), out_type(union_history.size());
  NumericVector out_from(union_history.size());
  CharacterVector input_reasons = unions["end_reason"];
  for (size_t u = 0; u < union_history.size(); ++u) {
    const UnionRecord& relationship = union_history[u];
    out_type[u] = relationship.cohab ? "cohab" : "marriage";
    out_from[u] = relationship.converted_from;
    out_uid[u] = relationship.id;
    out_wife[u] = people[relationship.wife].id;
    out_husband[u] = people[relationship.husband].id;
    out_start[u] = relationship.start;
    out_end[u] = relationship.end;
    if (NumericVector::is_na(relationship.end)) out_reason[u] = NA_STRING;
    else if (relationship.converted) out_reason[u] = "conversion";
    else if (relationship.dissolved) out_reason[u] = "dissolution";
    else if (u < static_cast<size_t>(unions.nrows())) {
      NumericVector input_ends = unions["end_time"];
      if (!NumericVector::is_na(input_ends[u])) out_reason[u] = input_reasons[u];
      else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
    } else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
  }
  DataFrame out_unions = DataFrame::create(_["union_id"] = out_uid,
    _["wife_id"] = out_wife, _["husband_id"] = out_husband,
    _["start_time"] = out_start, _["end_time"] = out_end, _["end_reason"] = out_reason);
  if (cohab_enabled || unions.containsElementNamed("union_type"))
    out_unions = DataFrame::create(_["union_id"] = out_uid, _["wife_id"] = out_wife,
      _["husband_id"] = out_husband, _["start_time"] = out_start, _["end_time"] = out_end,
      _["end_reason"] = out_reason, _["union_type"] = out_type,
      _["previous_cohab_id"] = out_from);
  return List::create(_["persons"] = DataFrame::create(
    _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
    _["death_time"] = out_death, _["mother_id"] = out_mother, _["father_id"] = out_father),
    _["transitions"] = DataFrame::create(_["person_id"] = wrap(tx_person), _["time"] = wrap(tx_time),
      _["group_origin"] = wrap(tx_origin), _["group_dest"] = wrap(tx_dest)),
    _["groups"] = out_group, _["unions"] = out_unions, _["ancestry"] = out_ancestry,
    _["births"] = birth_count, _["deaths"] = death_count,
    _["dissolutions"] = dissolution_count, _["conversions"] = conversion_count,
    _["weighted_infeasible"] = weighted_infeasible, _["draw_approximate"] = draw_approximate,
    _["unrealized"] = unrealized, _["cells_consistent"] = cells_consistent,
    _["hetfert"] = wrap(fertility_weights),
    _["diagnostic_log_weights"] = [&]() {
      // Rows are new unions in order; columns are group indices. Unions formed
      // without weights (unit weights) record zeros.
      NumericMatrix out(static_cast<int>(diagnostic_log_weights.size()), group_count);
      for (size_t u = 0; u < diagnostic_log_weights.size(); ++u)
        for (int g = 0; g < group_count; ++g)
          out(u, g) = diagnostic_log_weights[u].empty() ? 0 : diagnostic_log_weights[u][g];
      return out;
    }());
}

// Deterministic validation hook for contribution updates. Operations use zero-
// based person indices; no demographic event or RNG is involved.
// [[Rcpp::export]]
List inspect_intensity_tracker_cpp(NumericVector female_born, IntegerVector female_groups,
                                   DataFrame search_rates, DataFrame fertility_rates,
                                   NumericVector times, IntegerVector persons,
                                   CharacterVector actions, IntegerVector previous,
                                   NumericVector ready, int group_count,
                                   Nullable<NumericVector> search_origins = R_NilValue,
                                   Nullable<DataFrame> probability_rates = R_NilValue) {
  const int n = female_born.size(), ops = times.size();
  if (!ops || female_groups.size() != n || persons.size() != ops || actions.size() != ops ||
      previous.size() != ops || ready.size() != ops) stop("Invalid intensity-tracker inputs.");
  double lo = times[0], hi = times[ops - 1];
  for (double born : female_born) lo = std::min(lo, born);
  WeightCalibrator calibration(lo, hi, group_count, kWeightBinWidth, 0);
  FemaleIntensityTracker tracker(calibration);
  const bool stratified = search_rates.containsElementNamed("previous_union");
  const int categories = stratified ? 3 : 1;
  const std::array<std::string, 3> labels = {"none", "separated", "widowed"};
  std::vector<std::vector<Band>> searches(group_count * categories), fertility(group_count);
  for (int g = 0; g < group_count; ++g) {
    fertility[g] = read_bands(fertility_rates, "", "", g);
    for (int p = 0; p < categories; ++p)
      searches[g * categories + p] = read_bands(search_rates, "", "", g,
        stratified ? labels[p] : "", stratified && search_uses_duration(search_rates, labels[p]));
  }
  const bool split_types = probability_rates.isNotNull();
  std::vector<std::vector<Band>> probabilities(group_count * 3);
  if (split_types) {
    DataFrame table(probability_rates);
    const bool by_previous = table.containsElementNamed("previous_union");
    for (int g = 0; g < group_count; ++g) for (int p = 0; p < 3; ++p)
      probabilities[g * 3 + p] = read_bands(table, "", "", g, by_previous ? labels[p] : "");
  }
  NumericVector origins(ops, NA_REAL);
  if (search_origins.isNotNull()) origins = NumericVector(search_origins);
  if (origins.size() != ops) stop("Invalid search origins.");
  List snapshots(ops);
  IntegerVector pending(ops);
  std::vector<bool> active(n, false);
  for (int k = 0; k < ops; ++k) {
    if (!R_finite(times[k]) || (k && times[k] < times[k - 1])) stop("Unordered tracker times.");
    tracker.advance(times[k]);
    const std::string action = as<std::string>(actions[k]);
    const int person = persons[k];
    if (person < 0 || person >= n || previous[k] < 0 || previous[k] > 2)
      stop("Invalid intensity-tracker person or category.");
    if (action == "activate") {
      const int g = female_groups[person];
      if (!active[person]) calibration.add_woman(g, female_born[person], 1);
      active[person] = true;
      tracker.activate(person, g, female_born[person],
        &searches[g * categories + (stratified ? previous[k] : 0)], &fertility[g], ready[k], times[k], origins[k],
        split_types ? &probabilities[g * 3 + previous[k]] : nullptr);
    } else if (action == "remove") {
      tracker.remove(person);
      if (active[person]) calibration.add_woman(female_groups[person], female_born[person], -1);
      active[person] = false;
    }
    else if (action == "ready") tracker.set_ready(person, ready[k], times[k]);
    else if (action != "inspect") stop("Unknown intensity-tracker operation.");
    snapshots[k] = calibration.intensity_snapshot(split_types);
    pending[k] = tracker.pending();
  }
  return List::create(_["intensities"] = snapshots, _["pending"] = pending,
    _["bins"] = calibration.bins(), _["first"] = std::floor(lo / kWeightBinWidth));
}
