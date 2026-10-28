#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <functional>
#include <memory>
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
  // last_birth enforces spacing.
  double id, born, died, mother, father, last_birth;
  bool female;
  int group;
};
// person is a zero-based position in people, not the public person_id.
// A version invalidates queued events after a change in union status.
enum class EventType { death, birth, search, dissolution };
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
  std::vector<double> men_, women_, fertile_; // [group * bins + bin]
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
    fertile_ = men_;
  }
  int bins() const { return bins_; }
  // Register a utility table before counts are added; returns its index.
  int add_table(const LagTable* table) {
    tables_.push_back(table);
    sums_.emplace_back(static_cast<size_t>(groups_) * groups_ * bins_, 0.0);
    return static_cast<int>(tables_.size()) - 1;
  }
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
  void add_fertile(int group, double born, double delta) { fertile_[cell(group, born)] += delta; }

  // Expected male flow shares follow from each woman-bin's intensity and
  // her group-level choice shares b[g] S[w,g] / sum_h b[h] S[w,h]. Find log
  // weights so those shares meet the mirrored target (Newton on a convex
  // objective, warm-started). Empty or uninformative states return b = 1.
  WeightSolution solve(double now, const std::vector<const std::vector<Band>*>& hazards,
                       bool fertility, int table_id, const std::vector<double>& delta,
                       std::vector<double>& warm) {
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
          const double at_risk = fertility ? fertile_[i] : women_[i];
          if (at_risk > 0) {
            // Average hazard over the bin's whole age span, so narrow or
            // misaligned rate bands are not missed by a single centre age.
            const double oldest = now - static_cast<double>(first_ + k) * width_;
            const double hazard = rate_integral(*hazards[g], oldest - width_, oldest) / width_;
            if (hazard > 0) { rows.push_back({g, k, at_risk * hazard}); intensity[g] += at_risk * hazard; }
          }
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

// Internal entry point: calibrate weights for a supplied unpartnered state.
// [[Rcpp::export]]
List inspect_weighted_calibration_cpp(NumericVector female_born, IntegerVector female_groups,
                                      LogicalVector female_fertile, NumericVector male_born,
                                      IntegerVector male_groups, DataFrame rates, List matching,
                                      int group_count, double now, double cutoff, double width,
                                      bool fertility) {
  if (female_groups.size() != female_born.size() || female_fertile.size() != female_born.size() ||
      male_groups.size() != male_born.size() || group_count < 1)
    stop("Invalid weighted-calibration inputs.");
  double lo = now, hi = now;
  for (double born : female_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  for (double born : male_born) { lo = std::min(lo, born); hi = std::max(hi, born); }
  WeightCalibrator calibrator(lo, hi, group_count, width, cutoff);
  const LagTable table = build_lag_table(matching, group_count, calibrator.bins(), width);
  calibrator.add_table(&table);
  for (R_xlen_t i = 0; i < male_born.size(); ++i) calibrator.add_man(male_groups[i], male_born[i], 1);
  for (R_xlen_t i = 0; i < female_born.size(); ++i) {
    calibrator.add_woman(female_groups[i], female_born[i], 1);
    if (female_fertile[i]) calibrator.add_fertile(female_groups[i], female_born[i], 1);
  }
  std::vector<std::vector<Band>> bands(group_count);
  std::vector<const std::vector<Band>*> hazards(group_count);
  for (int g = 0; g < group_count; ++g) {
    bands[g] = read_bands(rates, "", "", g);
    hazards[g] = &bands[g];
  }
  std::vector<double> warm;
  const WeightSolution s = calibrator.solve(now, hazards, fertility, 0,
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
    incremental.add_fertile(female_groups[i], female_born[i], 1);
  }
  for (R_xlen_t i = 0; i < male_born.size(); ++i) incremental.add_man(male_groups[i], male_born[i], 1);
  std::vector<std::vector<Band>> bands(group_count);
  std::vector<const std::vector<Band>*> hazards(group_count);
  for (int g = 0; g < group_count; ++g) { bands[g] = read_bands(rates, "", "", g); hazards[g] = &bands[g]; }
  const std::vector<double> delta = as<std::vector<double>>(matching["delta"]);
  std::vector<double> warm;
  for (R_xlen_t i = 0; i < op_times.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    incremental.solve(op_times[i], hazards, false, 0, delta, warm);
    incremental.add_man(male_groups[op_men[i]], male_born[op_men[i]], op_deltas[i]);
  }
  const double end = op_times[op_times.size() - 1];
  WeightCalibrator rebuilt = incremental;
  rebuilt.refresh();
  const double difference = incremental.max_sum_difference(rebuilt, 0);
  std::vector<double> warm_a, warm_b;
  const WeightSolution a = incremental.solve(end, hazards, false, 0, delta, warm_a);
  const WeightSolution b = rebuilt.solve(end, hazards, false, 0, delta, warm_b);
  return List::create(_["max_sum_difference"] = difference,
                      _["incremental"] = wrap(a.log_weights), _["rebuilt"] = wrap(b.log_weights));
}

struct UnionRecord {
  double id, start, end;
  int wife, husband;
  bool wife_died = false;
  bool dissolved = false;
  bool dissolution_sampled = false;
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
                  DataFrame dissolution_rates) {
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
  MalePool pool(min_birth, start + years, group_count, kWeightBinWidth);
  for (int i = 0; i < static_cast<int>(people.size()); ++i)
    if (!people[i].female && NumericVector::is_na(people[i].died))
      pool.add(i, people[i].group, people[i].born);
  // Draw models for the union and (if separate) father matching models.
  // An empty population never draws; skip building custom-formula grids then.
  const DrawModel union_draw = people.empty() ? DrawModel() :
    make_draw_model(matching, group_count, start + years - min_birth);
  DrawModel father_draw_storage;
  const DrawModel* father_draw = &union_draw;
  if (!people.empty() && SEXP(father_matching) != SEXP(matching)) {
    father_draw_storage = make_draw_model(father_matching, group_count, start + years - min_birth);
    father_draw = &father_draw_storage;
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
  // weighted_choice state (D066). The calibrator exists only for grouped
  // populations; until it is created below, pool hooks record nothing.
  std::unique_ptr<WeightCalibrator> calibrator;
  std::vector<char> counted_woman, counted_fertile;
  std::multimap<double, int> spacing_entries;
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
  // weighted_choice counts living, unpartnered women; fertility-driven intensity
  // also requires the birth-spacing interval to have ended (queued re-entry).
  auto add_counted_woman = [&](int woman, double now) {
    if (!calibrator || !people[woman].female || !NumericVector::is_na(people[woman].died) ||
        current_partner[woman] >= 0 || counted_woman[woman]) return;
    counted_woman[woman] = 1;
    calibrator->add_woman(people[woman].group, people[woman].born, 1);
    const double ready = people[woman].last_birth + birth_interval;
    if (now >= ready) {
      counted_fertile[woman] = 1;
      calibrator->add_fertile(people[woman].group, people[woman].born, 1);
    } else spacing_entries.emplace(ready, woman);
  };
  auto remove_counted_woman = [&](int woman) {
    if (!calibrator || !counted_woman[woman]) return;
    counted_woman[woman] = 0;
    calibrator->add_woman(people[woman].group, people[woman].born, -1);
    if (counted_fertile[woman]) {
      counted_fertile[woman] = 0;
      calibrator->add_fertile(people[woman].group, people[woman].born, -1);
    }
  };
  auto counted_birth = [&](int woman, double now) {
    if (!calibrator || !counted_woman[woman]) return;
    if (counted_fertile[woman]) {
      counted_fertile[woman] = 0;
      calibrator->add_fertile(people[woman].group, people[woman].born, -1);
    }
    spacing_entries.emplace(now + birth_interval, woman);
  };
  auto release_spacing = [&](double now) {
    while (!spacing_entries.empty() && spacing_entries.begin()->first <= now) {
      const int woman = spacing_entries.begin()->second;
      spacing_entries.erase(spacing_entries.begin());
      if (counted_woman[woman] && !counted_fertile[woman] &&
          now >= people[woman].last_birth + birth_interval) {
        counted_fertile[woman] = 1;
        calibrator->add_fertile(people[woman].group, people[woman].born, 1);
      }
    }
  };
  std::vector<const std::vector<Band>*> search_hazards, fertility_hazards;
  LagTable union_table, father_table_storage;
  int union_table_id = 0, father_table_id = 0;
  // A NULL father model shares the general model and calibration table. When a
  // model prevents all exogamy, every woman sees only her own group's men, so
  // a group weight cannot change her choice: skip calibration.
  auto applies_weights = [](const List& model) {
    return !as<bool>(model["prevent_exogamy_all"]);
  };
  const bool union_weighted = applies_weights(matching);
  const bool father_weighted = applies_weights(father_matching);
  std::vector<double> union_warm, father_warm;
  double weighted_infeasible = 0;
  // Partner and father searches share eligibility and the exact draw; only
  // form_union changes partner state or removes the selected man.
  auto find_candidate = [&](int woman, double now, const List& choice_model,
                            const DrawModel& draw, bool father,
                            std::vector<double>* applied) -> int {
    if (pool.empty()) return -1;
    std::vector<double> log_weights(group_count, 0.0);
    if (calibrator && (father ? father_weighted : union_weighted)) {
      // Searches and on-demand unions use the union-formation calibration;
      // rate-mode nonmarital fathers have their own, driven by fertility.
      release_spacing(now);
      const bool fertility_driven = !rate_initiation || father;
      const WeightSolution weights = calibrator->solve(now,
        fertility_driven ? fertility_hazards : search_hazards, fertility_driven,
        father ? father_table_id : union_table_id,
        as<std::vector<double>>(choice_model["delta"]), father ? father_warm : union_warm);
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
    partner = find_candidate(woman, now, matching, union_draw, false,
                             diagnostics ? &applied : nullptr);
    if (partner >= 0) {
      if (next_union_id >= 9007199254740991.0 ||
          union_history.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Union history exceeds supported ID or row limits.");
      const int u = static_cast<int>(union_history.size());
      union_history.push_back({++next_union_id, now, NA_REAL, woman, partner});
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
    std::vector<Band> male, female, partnered, unpartnered, search;
  };
  const int inheritance_method = inheritance["method"];
  const double passing = inheritance["passing"];
  IntegerVector inheritance_ranks = inheritance["ranks"];
  IntegerVector ancestry_groups = inheritance["ancestry_groups"];
  Function inherit_custom = inheritance["custom"];
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
  // weighted_choice needs schedules and group indices, so it starts here.
  // One group means unit weights, so calibration is unnecessary.
  if ((union_weighted || father_weighted) && group_count > 1) {
    double min_birth = start;
    for (const Person& person : people) min_birth = std::min(min_birth, person.born);
    calibrator.reset(new WeightCalibrator(min_birth, start + years, group_count,
                                          kWeightBinWidth, min_male_age));
    for (int group = 0; group < group_count; ++group) {
      search_hazards.push_back(&schedules[group].search);
      fertility_hazards.push_back(&schedules[group].unpartnered);
    }
    if (union_weighted) {
      union_table = build_lag_table(matching, group_count, calibrator->bins(), kWeightBinWidth);
      union_table_id = father_table_id = calibrator->add_table(&union_table);
    }
    if (father_weighted && SEXP(father_matching) != SEXP(matching)) {
      father_table_storage = build_lag_table(father_matching, group_count, calibrator->bins(),
                                             kWeightBinWidth);
      father_table_id = calibrator->add_table(&father_table_storage);
    }
    counted_woman.assign(people.size(), 0);
    counted_fertile.assign(people.size(), 0);
    for (int man = 0; man < static_cast<int>(people.size()); ++man)
      if (pool.contains(man)) calibrator->add_man(people[man].group, people[man].born, 1);
    for (int woman = 0; woman < static_cast<int>(people.size()); ++woman)
      add_counted_woman(woman, start);
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
  const bool has_dissolution = dissolution_rates.nrows() > 0;
  std::vector<std::vector<Band>> dissolution_schedules;
  if (has_dissolution) {
    const int count = dissolution_rates.containsElementNamed(".hammel_group") ?
      group_count * group_count : 1;
    for (int pair = 0; pair < count; ++pair)
      dissolution_schedules.push_back(read_bands(dissolution_rates, "", "", pair));
  }
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
    double dissolution_time = std::numeric_limits<double>::infinity();
    if (has_dissolution && p.female && current_partner[i] >= 0) {
      UnionRecord& relationship = union_history[person_unions[i].back()];
      if (!relationship.dissolution_sampled) {
        const int pair = dissolution_schedules.size() == 1 ? 0 :
          p.group + group_count * people[current_partner[i]].group;
        relationship.dissolution_time = candidate(waiting_time(dissolution_schedules[pair],
          now - relationship.start, threshold()), now);
        relationship.dissolution_sampled = true;
      }
      // Retain the union clock through births, including a dissolution due now
      // after a tied partnered birth. It is not a new zero-wait random draw.
      dissolution_time = relationship.dissolution_time;
    }
    const double due = std::min(std::min(death_time, birth_time),
                                std::min(dissolution_time, search_time));
    // Exact same-person ties: death, birth, dissolution, then search (D078).
    const auto type = death_time <= due ? EventType::death :
      (birth_time <= due ? EventType::birth :
       (dissolution_time <= due ? EventType::dissolution : EventType::search));
    if (due < end) queue.emplace(due, Event{i, type, version});
  };
  // Initialize one next event for each living person. Known age changes are
  // already incorporated, so birthdays alone do not trigger resampling.
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    schedule(i, start);
  }
  double birth_count = 0, death_count = 0, dissolution_count = 0;
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
        add_counted_woman(relationship.wife, now);
        // Her old birth candidate is invalid once her husband has died.
        if ((rate_initiation || has_dissolution) && i == relationship.husband)
          schedule(relationship.wife, now);
      }
      if (!people[i].female) remove_male(i);
      else remove_counted_woman(i);
      ++death_count;
    } else if (event.type == EventType::dissolution) {
      const int partner = current_partner[i];
      if (partner < 0) stop("Dissolution requires an ongoing union.");
      UnionRecord& relationship = union_history[person_unions[i].back()];
      relationship.end = now;
      relationship.dissolved = true;
      current_partner[i] = current_partner[partner] = -1;
      add_male(partner);
      add_counted_woman(i, now);
      ++dissolution_count;
      // Only her state-dependent birth/search hazards change. The man's death
      // candidate remains valid, and the old episode no longer carries risk.
      schedule(i, now);
    } else if (event.type == EventType::search) {
      form_union(i, now);
      schedule(i, now);
    } else {
      double father = NA_REAL;
      int partner = current_partner[i];
      if (!rate_initiation && partner < 0) partner = form_union(i, now);
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
      counted_birth(i, now);
      // Default ffab is 0.4886 - its the signature of Berkeley training!
      bool female = R::runif(0, 1) < ffab;
      int child = static_cast<int>(people.size());
      people.push_back({++next_id, now, NA_REAL, mother, father,
                        -std::numeric_limits<double>::infinity(), female, people[i].group});
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
      if (!female) {
        pool.add(child, child_group, now);
        if (calibrator) calibrator->add_man(child_group, now, 1);
      }
      person_unions.emplace_back();
      current_partner.push_back(-1);
      if (calibrator) {
        counted_woman.push_back(0);
        counted_fertile.push_back(0);
        add_counted_woman(child, now);
      }
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
  return List::create(_["persons"] = DataFrame::create(
    _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
    _["death_time"] = out_death, _["mother_id"] = out_mother, _["father_id"] = out_father),
    _["groups"] = out_group, _["unions"] = out_unions, _["ancestry"] = out_ancestry,
    _["births"] = birth_count, _["deaths"] = death_count,
    _["dissolutions"] = dissolution_count,
    _["weighted_infeasible"] = weighted_infeasible, _["draw_approximate"] = draw_approximate,
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
