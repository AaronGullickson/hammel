#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <functional>
#include <memory>
#include <set>
#include <tuple>
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
// Shared by the calibration tracker and fertility/mortality cells: entries are
// rate changes (age bands, gaps, spacing release), not demographic events, and
// consume no randomness. Ties are ordered by person index.
//
// A bucket queue: a ring of buckets of width 1/kBuckets years spanning
// kSpan years, plus an overflow heap for later entries. The span covers the
// usual band widths (input bands of five or more years for rate processes), so
// few boundaries need the heap. A bucket is sorted by
// (time, person) when reached, so pops follow exactly (time, person) order.
// Entries are lazily deleted: one is current only while it matches next_[person].
class BoundaryQueue {
  static constexpr int kBuckets = 365, kSpan = 16, kRing = kBuckets * kSpan;
  using Entry = std::pair<double, int>;
  struct Later {
    bool operator()(const Entry& a, const Entry& b) const { return a > b; }
  };
  std::vector<std::vector<Entry>> ring_ = std::vector<std::vector<Entry>>(kRing);
  std::vector<Entry> overflow_;          // min-heap of entries beyond the ring
  std::vector<double> next_;
  size_t live_ = 0, in_ring_ = 0, pos_ = 0;
  long long cur_ = 0;                    // absolute index of the bucket being drained
  bool started_ = false, sorted_ = false;
  static long long bucket(double t) { return static_cast<long long>(std::floor(t * kBuckets)); }
  std::vector<Entry>& slot(long long b) {
    return ring_[static_cast<size_t>(((b % kRing) + kRing) % kRing)];
  }
  void push(const Entry& e) {
    long long b = bucket(e.first);
    if (!started_) { cur_ = b; started_ = true; }
    if (b < cur_) b = cur_;              // at or after the current time; roundoff guard
    if (b - cur_ >= kRing) {
      overflow_.push_back(e);
      std::push_heap(overflow_.begin(), overflow_.end(), Later());
      return;
    }
    auto& v = slot(b);
    if (b == cur_ && sorted_)
      v.insert(std::upper_bound(v.begin() + static_cast<long>(pos_), v.end(), e), e);
    else v.push_back(e);
    ++in_ring_;
  }
  // Position the cursor on the earliest current entry; false if none remain.
  bool settle() {
    if (live_ == 0) return false;
    for (;;) {
      if (!sorted_ && in_ring_ == 0 && !overflow_.empty()) // jump over empty buckets
        cur_ = std::max(cur_, bucket(overflow_.front().first));
      auto& v = slot(cur_);
      if (!sorted_) {
        while (!overflow_.empty() && bucket(overflow_.front().first) <= cur_) {
          v.push_back(overflow_.front());
          ++in_ring_;
          std::pop_heap(overflow_.begin(), overflow_.end(), Later());
          overflow_.pop_back();
        }
        std::sort(v.begin(), v.end());
        sorted_ = true; pos_ = 0;
      }
      while (pos_ < v.size()) {
        const Entry& e = v[pos_];
        if (e.first == next_[e.second]) return true;
        ++pos_; --in_ring_;
      }
      v.clear(); sorted_ = false; ++cur_;
    }
  }
 public:
  void set(int person, double time) {
    if (person >= static_cast<int>(next_.size())) next_.resize(person + 1, R_PosInf);
    if (next_[person] == time) return;
    if (R_finite(next_[person])) --live_;
    next_[person] = time;
    if (R_finite(time)) { ++live_; push({time, person}); }
  }
  void clear(int person) {
    if (person < static_cast<int>(next_.size())) set(person, R_PosInf);
  }
  double top() {
    if (!settle()) return R_PosInf;
    return slot(cur_)[pos_].first;
  }
  // Remove and return the earliest entry; the caller refreshes that person.
  Entry pop() {
    settle();
    const Entry e = slot(cur_)[pos_++];
    --in_ring_;
    next_[e.second] = R_PosInf;
    --live_;
    return e;
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
    if (on_probability_) on_probability_(person, p);
  }
  // Weighted search (step 2) keeps search intensities elsewhere and needs each
  // counted woman's current cohabitation probability.
  std::function<void(int, double)> on_probability_;
 public:
  explicit FemaleIntensityTracker(WeightCalibrator& calibration) : calibration_(calibration) {}
  void set_probability_hook(std::function<void(int, double)> f) { on_probability_ = std::move(f); }
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
// keeps its intensity: when it is drawn, no child is born and one unrealized birth is recorded.
// Clocks keep an Exp(1) threshold consumed at the current intensity, so count
// changes recompute firing times without drawing random numbers.
// Upper-open fertility bands stay one cell when birth_interval * rate is below
// this share (a numerical tolerance, not a demographic parameter).
constexpr double kOpenBandSplitTolerance = 0.001;

// Sum tree over a process's cells: leaf c holds that cell's current event
// intensity. Sums are recomputed from children along each updated path, so
// they do not drift. Capacity grows geometrically as cells are added.
class IntensityTree {
  size_t capacity_ = 0;
  std::vector<double> tree_;
 public:
  void set(size_t cell, double value) {
    if (cell >= capacity_) {
      size_t next = capacity_ ? capacity_ : 1;
      while (next <= cell) next *= 2;
      std::vector<double> grown(2 * next, 0);
      for (size_t i = 0; i < capacity_; ++i) grown[next + i] = tree_[capacity_ + i];
      for (size_t i = next - 1; i > 0; --i) grown[i] = grown[2 * i] + grown[2 * i + 1];
      capacity_ = next;
      tree_.swap(grown);
    }
    size_t node = capacity_ + cell;
    tree_[node] = value;
    while (node > 1) {
      node /= 2;
      tree_[node] = tree_[2 * node] + tree_[2 * node + 1];
    }
  }
  double total() const { return tree_.empty() ? 0 : tree_[1]; }
  double value(size_t cell) const { return cell < capacity_ ? tree_[capacity_ + cell] : 0; }
  // The cell holding position uniform * total; only cells with positive intensity.
  size_t choose(double uniform) const {
    double value = std::min(uniform * total(), std::nextafter(total(), 0.0));
    size_t node = 1;
    while (node < capacity_) {
      const double left = tree_[2 * node], right = tree_[2 * node + 1];
      // Roundoff guard: never descend into a branch with no intensity.
      if ((value < left && left > 0) || !(right > 0)) node *= 2;
      else { value = std::max(0.0, value - left); node = 2 * node + 1; }
    }
    return node - capacity_;
  }
};

// One event clock for a whole process. Independent cell streams sum to one
// Poisson stream at the total intensity; the cell of each event is then chosen
// in proportion to its intensity. The Exp(1) threshold is spent at the current
// total, so intensity changes recompute the due time without random draws.
// Call spend(now) before every intensity change.
class ProcessClock {
  double threshold_ = 0, last_ = 0;
  bool has_threshold_ = false;
 public:
  void start(double now) { last_ = now; }
  void spend(const IntensityTree& tree, double now) {
    if (has_threshold_ && now > last_)
      threshold_ = std::max(0.0, threshold_ - tree.total() * (now - last_));
    last_ = now;
  }
  // After intensity changes: ensure a threshold exists while events can occur.
  void arm(const IntensityTree& tree, const std::function<double()>& draw) {
    if (!has_threshold_ && tree.total() > 0) { threshold_ = draw(); has_threshold_ = true; }
  }
  double due(const IntensityTree& tree) const {
    return has_threshold_ && tree.total() > 0 ? last_ + threshold_ / tree.total() : R_PosInf;
  }
  // The clock fired at now: its threshold is used up.
  void fired(double now) { has_threshold_ = false; last_ = now; }
};

// Choose a member from log weights with one uniform draw. Equal log weights
// use the unweighted rule floor(uniform * n), so constant weights reproduce
// unweighted runs exactly (milestone 8 step 2).
size_t pick_log_weights(const std::vector<double>& eta, double uniform) {
  const size_t n = eta.size();
  bool equal = true;
  double top = eta[0];
  for (double e : eta) { if (e != eta[0]) equal = false; top = std::max(top, e); }
  if (equal) return std::min(n - 1, static_cast<size_t>(std::floor(uniform * static_cast<double>(n))));
  std::vector<double> cumulative(n);
  double total = 0;
  for (size_t k = 0; k < n; ++k) { total += std::exp(eta[k] - top); cumulative[k] = total; }
  const double target = uniform * total;
  for (size_t k = 0; k < n; ++k) if (target < cumulative[k]) return k;
  return n - 1;
}

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
  void update(size_t slot, double value) { set(slot, value); }
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
    int group = 0;             // her group (cells pool groups when rates do)
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
    double lambda = 0;
    std::map<size_t, std::pair<int, FertilityWeightSum>> counted; // count and weight per unpartnered calibration bin
    std::map<int, int> by_group;               // members by group, when groups share cells
                                               // (attributes unrealized births to a group)
  };
  struct Member {
    int cell = -1, slot = -1, group = 0;
    bool eligible = false, counted = false;
    size_t bin = 0;
    double weight = 1;
  };
  std::vector<Source> sources_;
  std::vector<Stratum> strata_;
  std::vector<Cell> cells_;
  std::vector<Member> members_;
  std::vector<State> states_;   // engine state at each woman's last refresh
  // Weight models (step 2). Without matching calibration, log weights are
  // computed for the chosen cell's eligible women at selection. With
  // calibration, refresh_eta_ folds them into each woman's stored weight.
  std::function<void(int, const std::vector<int>&, std::vector<double>&)> select_eta_;
  std::function<double(int, int)> refresh_eta_;
  // lookup_[group * types + type][actual parity] -> stratum (-1 until created);
  // without parity-specific schedules, index 0 serves every parity.
  int types_;
  bool split_parity_, weighted_;
  bool pooled_groups_ = false;   // rates shared by groups: one set of cells
  double max_age_, interval_;
  std::vector<std::vector<int>> lookup_;
  BoundaryQueue boundaries_;
  // One birth clock for all cells; every cell contributes rate * members.
  IntensityTree intensity_;
  ProcessClock clock_;
  std::function<State(int)> state_;
  std::function<double()> draw_threshold_;
  std::vector<double> eta_;     // scratch log weights for selection
  const WeightCalibrator* calibrator_;
  double start_, end_, period_origin_;
  // Called with (cell, group) when a cell without eligible women is drawn.
  std::function<void(int, int)> on_unrealized_;
  unsigned long long updates_ = 0;

  // Spend the clock for the interval since the last change. Must run before
  // membership or eligibility changes.
  void settle(int, double now) { clock_.spend(intensity_, now); }
  // Recompute the cell's intensity from current counts (after settle).
  void touch(int c, double now) {
    Cell& x = cells_[c];
    x.lambda = x.rate * x.members;
    clock_.spend(intensity_, now);
    intensity_.set(static_cast<size_t>(c), x.lambda);
    clock_.arm(intensity_, draw_threshold_);
  }
  void leave(int woman, int& changed) {
    Member& m = members_[woman];
    if (m.cell < 0) return;
    Cell& x = cells_[m.cell];
    --x.members;
    if (pooled_groups_ && --x.by_group[m.group] == 0) x.by_group.erase(m.group);
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
      cells_.push_back(x);
    }
    strata_.push_back(st);
    return k;
  }

 public:
  FertilityCells(std::vector<Source> sources, int groups, int types, bool split_parity,
                 double max_age, double interval, bool weighted,
                 std::function<State(int)> state, std::function<double()> draw_threshold,
                 const WeightCalibrator* calibrator, double start, double end, double period_origin)
      : sources_(std::move(sources)), types_(types), split_parity_(split_parity), weighted_(weighted),
        max_age_(max_age), interval_(interval), lookup_(static_cast<size_t>(groups) * types),
        state_(std::move(state)), draw_threshold_(std::move(draw_threshold)),
        calibrator_(calibrator), start_(start), end_(end), period_origin_(period_origin) {
    clock_.start(start);
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
    if (woman >= static_cast<int>(members_.size())) {
      members_.resize(woman + 1);
      states_.resize(woman + 1);
    }
    states_[woman] = state_(woman);
    place(woman, now);
  }
  // Place a woman from her cached state (set by the last refresh). Only time
  // changes between refreshes, so boundaries reuse the cached state.
  void place(int woman, double now) {
    const State& s = states_[woman];
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
    const size_t bin = counted ? calibrator_->bin_index(s.group, s.born) : 0;
    double weight = s.weight;
    if (refresh_eta_ && cell >= 0) {
      weight *= std::exp(refresh_eta_(woman, cell));
      if (!(weight > 0) || !R_finite(weight)) stop("A fertility weight is not finite and positive.");
    }
    Member& m = members_[woman];
    if (m.cell != cell || m.eligible != eligible || m.counted != counted || (counted && m.bin != bin) ||
        (cell >= 0 && weight != m.weight)) {
      // Settle elapsed intensity under the old state before changing it.
      if (m.cell >= 0) settle(m.cell, now);
      if (cell >= 0 && cell != m.cell) settle(cell, now);
      int old_cell = -1;
      leave(woman, old_cell);
      if (cell >= 0) {
        Cell& x = cells_[cell];
        ++x.members;
        if (pooled_groups_) ++x.by_group[s.group];
        m.cell = cell; m.eligible = eligible; m.counted = counted; m.bin = bin; m.weight = weight;
        m.group = s.group;
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
    place(entry.second, entry.first);
    if (++updates_ % 1024 == 0) checkUserInterrupt();
  }
  double next_clock() const {
    const double due = clock_.due(intensity_);
    return due < end_ ? due : R_PosInf;
  }
  // The birth clock fires: choose a cell by intensity, then a mother in it.
  int fire_next(double now) {
    clock_.spend(intensity_, now);
    const int c = static_cast<int>(intensity_.choose(R::runif(0, 1)));
    Cell& x = cells_[c];
    const double uniform = R::runif(0, 1);
    if (x.eligible.empty()) {
      // No eligible woman: an unrealized birth. In a cell shared by groups the
      // uniform draw picks the group in proportion to its members.
      int group = strata_[x.stratum].group;
      if (pooled_groups_) {
        double target = uniform * static_cast<double>(x.members);
        for (const auto& g : x.by_group) {
          group = g.first;
          if ((target -= g.second) < 0) break;
        }
      }
      clock_.fired(now);
      clock_.arm(intensity_, draw_threshold_);
      if (on_unrealized_) on_unrealized_(c, group);
      return -1;
    }
    const double n = static_cast<double>(x.eligible.size());
    size_t pick = weighted_ ? x.weights.choose(uniform) :
      std::min(x.eligible.size() - 1, static_cast<size_t>(std::floor(uniform * n)));
    if (select_eta_) {
      // Equal formula weights keep the pick above (hetfert or uniform).
      select_eta_(c, x.eligible, eta_);
      bool equal = true;
      for (double e : eta_) if (e != eta_[0]) { equal = false; break; }
      if (!equal) {
        if (weighted_) for (size_t k = 0; k < eta_.size(); ++k)
          eta_[k] += std::log(members_[x.eligible[k]].weight);
        pick = pick_log_weights(eta_, uniform);
      }
    }
    const int woman = x.eligible[pick];
    clock_.fired(now);
    clock_.arm(intensity_, draw_threshold_);
    return woman;
  }
  void set_select_eta(std::function<void(int, const std::vector<int>&, std::vector<double>&)> f) {
    select_eta_ = std::move(f);
  }
  void set_refresh_eta(std::function<double(int, int)> f) { refresh_eta_ = std::move(f); }
  void set_pooled_groups(bool pooled) { pooled_groups_ = pooled; }
  // A cell's group, union-type index, parity range and age band.
  struct CellInfo { int group, type; double parity_start, parity_end, age_start, age_end; };
  CellInfo cell_info(int c) const {
    const Cell& x = cells_[c];
    const Stratum& st = strata_[x.stratum];
    return {st.group, st.type, st.parity_start, st.parity_end, x.age_start, x.age_end};
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
  void set_unrealized_hook(std::function<void(int, int)> f) { on_unrealized_ = std::move(f); }
  // Calendar rate changes: update rates in place when the new sources keep
  // every stratum's cells (same sources and, after age splitting, the same
  // band edges). Returns false, changing nothing, otherwise.
  bool update_rates(const std::vector<Source>& sources, double now) {
    if (sources.size() != sources_.size()) return false;
    for (size_t i = 0; i < sources.size(); ++i)
      if (sources[i].group != sources_[i].group || sources[i].type != sources_[i].type ||
          sources[i].parity_start != sources_[i].parity_start || sources[i].parity_end != sources_[i].parity_end)
        return false;
    std::vector<std::vector<Band>> staged(strata_.size());
    for (size_t k = 0; k < strata_.size(); ++k) {
      const Stratum& st = strata_[k];
      const Source* source = nullptr;
      for (const Source& x : sources)
        if (x.group == st.group && x.type == st.type && x.parity_start <= st.parity_start &&
            st.parity_start < x.parity_end) { source = &x; break; }
      if (!source) return false;
      staged[k] = split_ages(source->bands, max_age_, interval_);
      if (staged[k].size() != st.bands.size()) return false;
      for (size_t i = 0; i < staged[k].size(); ++i)
        if (staged[k][i].start != st.bands[i].start || staged[k][i].end != st.bands[i].end) return false;
    }
    clock_.spend(intensity_, now);
    for (size_t k = 0; k < strata_.size(); ++k) {
      Stratum& st = strata_[k];
      for (size_t i = 0; i < staged[k].size(); ++i) {
        const int c = st.first_cell + static_cast<int>(i);
        st.bands[i].rate = staged[k][i].rate;
        Cell& x = cells_[c];
        x.rate = staged[k][i].rate;
        x.lambda = x.rate * x.members;
        intensity_.set(static_cast<size_t>(c), x.lambda);
      }
    }
    clock_.arm(intensity_, draw_threshold_);
    sources_ = sources;
    return true;
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
    for (size_t c = 0; c < cells_.size(); ++c) {
      if (members[c] != cells_[c].members || eligible[c] != static_cast<int>(cells_[c].eligible.size()))
        return false;
      if (intensity_.value(c) != cells_[c].lambda) return false;
    }
    if (weighted_) for (size_t c = 0; c < cells_.size(); ++c)
      if (std::abs(weights[c] - cells_[c].weights.sum()) > 1e-9 * std::max(weights[c], 1.0)) return false;
    return true;
  }
};

// ---- Cell-scheduled rate processes (milestone 8 step 1) --------------------
// One demographic rate process (mortality, union search, dissolution,
// conversion, or one group-transition slot). Each at-risk unit (a person, or a
// union identified by its history index) belongs to one cell: one input band
// of its schedule, measured from the unit's own clock origin (birth, union
// start, last union end, or group entry). Cells are the input bands: with no
// eligibility rules every member's hazard is the cell rate, so splitting a band
// cannot change results, only add boundary updates. Zero-rate bands are gaps.
// A cell generates events at rate * members; one process clock fires at the
// total intensity, the cell is chosen in proportion to its intensity, and the
// unit uniformly within it. No eligibility rules or weights.
class RateCells {
  struct Cell {
    int stratum = 0;
    double rate = 0, start = 0, end = 0;   // band edges on the unit's clock
    std::vector<int> members;   // swap-remove list of units
    double weight_bound = 0;     // upper bound on members' stored weights
  };
  struct Member { int cell = -1, slot = -1, stratum = -1; double origin = 0, weight = 1; };
  bool stored_weights_ = false;
  // strata_[k] = split bands of schedule k; first_[k] = its first cell index.
  std::vector<std::vector<Band>> strata_;
  std::vector<int> first_;
  std::vector<Cell> cells_;
  std::vector<Member> members_;
  BoundaryQueue boundaries_;
  // One event clock for all cells (see ProcessClock).
  IntensityTree intensity_;
  ProcessClock clock_;
  std::function<double()> draw_threshold_;
  // Weighted selection (step 2): returns the index of the chosen member of a
  // cell's member list, given the cell and one uniform draw.
  std::function<size_t(int, const std::vector<int>&, double)> chooser_;
  // Called after a unit changes cell (old, new; -1 for none).
  std::function<void(int, int, int)> on_move_;
  double end_;
  unsigned long long updates_ = 0;
  void update(int c, double now) {
    clock_.spend(intensity_, now);
    intensity_.set(static_cast<size_t>(c),
                   cells_[c].rate * static_cast<double>(cells_[c].members.size()));
    clock_.arm(intensity_, draw_threshold_);
  }

 public:
  RateCells(const std::vector<std::vector<Band>>& schedules,
                 std::function<double()> draw_threshold, double start, double end)
      : draw_threshold_(std::move(draw_threshold)), end_(end) {
    for (const auto& bands : schedules) {
      std::vector<Band> cells;
      for (const Band& band : bands) if (band.rate > 0) cells.push_back(band);
      strata_.push_back(std::move(cells));
      first_.push_back(static_cast<int>(cells_.size()));
      for (const Band& band : strata_.back()) {
        Cell x;
        x.stratum = static_cast<int>(strata_.size()) - 1;
        x.rate = band.rate; x.start = band.start; x.end = band.end;
        cells_.push_back(std::move(x));
      }
    }
    if (!cells_.empty()) intensity_.set(cells_.size() - 1, 0);
    clock_.start(start);
  }
  // Place a unit in the cell of schedule `stratum` (-1: not at risk) that
  // holds its clock value now - origin, and queue its next band boundary.
  void refresh(int person, int stratum, double origin, double now) {
    if (person >= static_cast<int>(members_.size())) members_.resize(person + 1);
    if (stratum >= 0 && !R_finite(origin)) stop("A rate process requires a finite clock origin.");
    // Fast path: most refreshes follow events that leave this process's state
    // unchanged. Same stratum and origin, still inside the current cell's band:
    // the cell and the queued boundary are already correct.
    Member& m = members_[person];
    if (stratum == m.stratum && (stratum < 0 ? m.cell < 0 : origin == m.origin && m.cell >= 0 &&
        origin + cells_[m.cell].start <= now && now < origin + cells_[m.cell].end)) return;
    int cell = -1;
    double next = R_PosInf;
    if (stratum >= 0) {
      const auto& bands = strata_[stratum];
      size_t lo = 0, hi = bands.size();
      while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (origin + bands[mid].start <= now) lo = mid + 1; else hi = mid;
      }
      if (lo > 0 && now < origin + bands[lo - 1].end) {
        cell = first_[stratum] + static_cast<int>(lo - 1);
        next = origin + bands[lo - 1].end;
      } else if (lo < bands.size()) {
        next = origin + bands[lo].start;
      }
    }
    m.stratum = stratum; m.origin = origin;
    if (m.cell != cell) {
      const int old = m.cell;
      if (old >= 0) {
        auto& list = cells_[old].members;
        const int last = list.back();
        list[m.slot] = last;
        members_[last].slot = m.slot;
        list.pop_back();
        if (list.empty()) cells_[old].weight_bound = 0;
      }
      m.cell = cell;
      if (cell >= 0) {
        m.slot = static_cast<int>(cells_[cell].members.size());
        cells_[cell].members.push_back(person);
        if (stored_weights_) cells_[cell].weight_bound = std::max(cells_[cell].weight_bound, m.weight);
      } else m.slot = -1;
      if (old >= 0) update(old, now);
      if (cell >= 0) update(cell, now);
      if (on_move_) on_move_(person, old, cell);
    }
    boundaries_.set(person, next);
  }
  double next_boundary() { return boundaries_.top(); }
  // Move the unit with the earliest boundary to its next cell, using the
  // schedule and clock origin recorded at its last refresh.
  void advance_boundary() {
    if (++updates_ % 1024 == 0) checkUserInterrupt();
    const auto entry = boundaries_.pop();
    const Member& m = members_[entry.second];
    refresh(entry.second, m.stratum, m.origin, entry.first);
  }
  double next_clock() const {
    const double due = clock_.due(intensity_);
    return due < end_ ? due : R_PosInf;
  }
  // The clock fires: choose a cell by intensity, then a member uniformly, and
  // re-arm. The caller updates the unit's state and refreshes its membership.
  int fire_next(double now) {
    clock_.spend(intensity_, now);
    const int c = static_cast<int>(intensity_.choose(R::runif(0, 1)));
    const Cell& x = cells_[c];
    if (x.members.empty()) stop("A rate cell fired without members.");
    const double n = static_cast<double>(x.members.size());
    const double uniform = R::runif(0, 1);
    size_t pick = std::min(x.members.size() - 1, static_cast<size_t>(std::floor(uniform * n)));
    if (stored_weights_) {
      // Stored weights: rejection sampling against the cell's weight bound,
      // exact for any bound at or above the largest weight. A member at the
      // bound is accepted without a further draw, so equal weights keep the
      // unweighted pick and its random-number use.
      for (;;) {
        const double w = members_[x.members[pick]].weight;
        if (w >= x.weight_bound || R::runif(0, 1) * x.weight_bound < w) break;
        pick = std::min(x.members.size() - 1, static_cast<size_t>(std::floor(R::runif(0, 1) * n)));
      }
    } else if (chooser_) pick = chooser_(c, x.members, uniform);
    clock_.fired(now);
    clock_.arm(intensity_, draw_threshold_);
    return x.members[pick];
  }
  void set_chooser(std::function<size_t(int, const std::vector<int>&, double)> f) {
    chooser_ = std::move(f);
  }
  void set_move_hook(std::function<void(int, int, int)> f) { on_move_ = std::move(f); }
  // Calendar rate changes: true if new schedules keep this process's cells
  // (same strata and nonzero band edges), so rates can change in place.
  bool same_layout(const std::vector<std::vector<Band>>& schedules) const {
    if (schedules.size() != strata_.size()) return false;
    for (size_t k = 0; k < schedules.size(); ++k) {
      size_t i = 0;
      for (const Band& band : schedules[k]) {
        if (!(band.rate > 0)) continue;
        if (i >= strata_[k].size() || band.start != strata_[k][i].start || band.end != strata_[k][i].end)
          return false;
        ++i;
      }
      if (i != strata_[k].size()) return false;
    }
    return true;
  }
  // Change every cell's rate at now (after same_layout). The clock spends its
  // threshold under the old rates and continues; nobody changes cells.
  void set_rates(const std::vector<std::vector<Band>>& schedules, double now) {
    clock_.spend(intensity_, now);
    for (size_t k = 0; k < schedules.size(); ++k) {
      size_t i = 0;
      for (const Band& band : schedules[k]) {
        if (!(band.rate > 0)) continue;
        const int c = first_[k] + static_cast<int>(i);
        strata_[k][i].rate = band.rate;
        cells_[c].rate = band.rate;
        intensity_.set(static_cast<size_t>(c), band.rate * static_cast<double>(cells_[c].members.size()));
        ++i;
      }
    }
    clock_.arm(intensity_, draw_threshold_);
  }
  // Keep a weight per member; events choose members in proportion to it.
  void enable_stored_weights() { stored_weights_ = true; }
  void set_weight(int unit, double weight) {
    if (unit >= static_cast<int>(members_.size())) members_.resize(unit + 1);
    Member& m = members_[unit];
    if (m.weight == weight) return;
    m.weight = weight;
    if (m.cell >= 0) cells_[m.cell].weight_bound = std::max(cells_[m.cell].weight_bound, weight);
  }
  int cell_of(int unit) const {
    return unit < static_cast<int>(members_.size()) ? members_[unit].cell : -1;
  }
  double cell_rate(int c) const { return cells_[c].rate; }
  // A cell's stratum and band edges on its clock.
  int cell_stratum(int c) const { return cells_[c].stratum; }
  double cell_start(int c) const { return cells_[c].start; }
  double cell_end(int c) const { return cells_[c].end; }
  bool consistent() const {
    std::vector<int> counts(cells_.size(), 0);
    for (size_t p = 0; p < members_.size(); ++p) {
      const Member& m = members_[p];
      if (m.cell < 0) continue;
      ++counts[m.cell];
      if (cells_[m.cell].members[m.slot] != static_cast<int>(p)) return false;
    }
    if (stored_weights_) {
      // Stored weights: each cell's bound covers its members' weights.
      for (size_t c = 0; c < cells_.size(); ++c)
        for (int unit : cells_[c].members)
          if (members_[unit].weight > cells_[c].weight_bound) return false;
    }
    for (size_t c = 0; c < cells_.size(); ++c) {
      if (counts[c] != static_cast<int>(cells_[c].members.size())) return false;
      if (intensity_.value(c) != cells_[c].rate * static_cast<double>(counts[c])) return false;
    }
    return true;
  }
};

// ---- Weighted search intensities for matching calibration (step 2) --------
// With search weights, a woman's search hazard is rate_c * N_c * w_i / W_c, so
// her calibration contribution depends on her whole cell. Each search cell
// keeps, per calibration bin, the member count and the sums of w (1 - p) and
// w p (p: her cohabitation probability). A bin's marriage and cohabitation
// intensities are sum over cells of rate_c * N_c / W_c times those sums.
// Cells whose contents changed are recomputed before the next search solve.
class SearchCalibration {
  struct Woman { int cell = -1; size_t bin = 0; double weight = 1; bool active = false; };
  // Per calibration bin within a cell: members, sums, and the last applied
  // contribution (kept until applied as zero, so stale entries are removed).
  struct BinSums { int count = 0; double marriage = 0, cohab = 0; std::array<double, 2> applied = {0, 0}; };
  struct CellSums { int count = 0; double weight = 0; std::map<size_t, BinSums> bins; bool dirty = false; };
  std::vector<Woman> women_;
  std::vector<double> probability_;
  std::vector<CellSums> cells_;
  std::vector<int> dirty_;
  std::vector<std::array<double, 2>> totals_;
  std::vector<std::array<int, 2>> positive_;
  const RateCells* search_ = nullptr;
  void grow(int woman) {
    if (woman >= static_cast<int>(women_.size())) {
      women_.resize(woman + 1);
      probability_.resize(woman + 1, 0.0);
    }
  }
  CellSums& cell(int c) {
    if (c >= static_cast<int>(cells_.size())) cells_.resize(c + 1);
    return cells_[c];
  }
  void touch(int c) {
    CellSums& x = cell(c);
    if (!x.dirty) { x.dirty = true; dirty_.push_back(c); }
  }
  void remove(int woman) {
    Woman& x = women_[woman];
    if (!x.active) return;
    CellSums& c = cell(x.cell);
    const double p = probability_[woman];
    BinSums& b = c.bins[x.bin];
    b.marriage -= x.weight * (1 - p);
    b.cohab -= x.weight * p;
    if (--b.count == 0) b.marriage = b.cohab = 0;
    c.weight -= x.weight;
    if (--c.count == 0) c.weight = 0;
    touch(x.cell);
    x.active = false;
  }
  void add(int woman) {
    Woman& x = women_[woman];
    CellSums& c = cell(x.cell);
    const double p = probability_[woman];
    BinSums& b = c.bins[x.bin];
    b.marriage += x.weight * (1 - p);
    b.cohab += x.weight * p;
    ++b.count;
    c.weight += x.weight;
    ++c.count;
    touch(x.cell);
    x.active = true;
  }

 public:
  void reset(const RateCells* search, size_t bins) {
    search_ = search;
    for (Woman& x : women_) x.active = false;
    cells_.clear();
    dirty_.clear();
    totals_.assign(bins, {0, 0});
    positive_.assign(bins, {0, 0});
  }
  // Place a woman in a search cell (-1: not searching) with her bin and weight.
  void update(int woman, int c, size_t bin, double weight) {
    grow(woman);
    const Woman& x0 = women_[woman];
    if (x0.active ? (c == x0.cell && bin == x0.bin && weight == x0.weight) : c < 0) return;
    remove(woman);
    if (c < 0) return;
    Woman& x = women_[woman];
    x.cell = c; x.bin = bin; x.weight = weight;
    add(woman);
  }
  void set_probability(int woman, double p) {
    grow(woman);
    if (probability_[woman] == p) return;
    const bool active = women_[woman].active;
    remove(woman);
    probability_[woman] = p;
    if (active) add(woman);
  }
  double weight(int woman) const { return women_[woman].weight; }
  // Search rates changed in place: every cell's contribution is recomputed.
  void mark_all_dirty() { for (size_t c = 0; c < cells_.size(); ++c) touch(static_cast<int>(c)); }
  size_t bin(int woman) const { return women_[woman].bin; }
  bool active(int woman) const { return woman < static_cast<int>(women_.size()) && women_[woman].active; }
  // Recompute changed cells and return marriage and cohabitation totals by bin.
  void apply(std::vector<double>& marriage, std::vector<double>& cohab) {
    for (int c : dirty_) {
      CellSums& x = cells_[c];
      x.dirty = false;
      const double scale = x.count > 0 && x.weight > 0 ?
        search_->cell_rate(c) * static_cast<double>(x.count) / x.weight : 0.0;
      for (auto it = x.bins.begin(); it != x.bins.end();) {
        BinSums& e = it->second;
        const size_t b = it->first;
        const std::array<double, 2> value = {e.count ? scale * std::max(0.0, e.marriage) : 0.0,
                                             e.count ? scale * std::max(0.0, e.cohab) : 0.0};
        for (int k = 0; k < 2; ++k) {
          positive_[b][k] += (value[k] > 0) - (e.applied[k] > 0);
          totals_[b][k] += value[k] - e.applied[k];
          // Exact zero when no cell contributes: no phantom roundoff intensity.
          if (positive_[b][k] == 0) totals_[b][k] = 0;
        }
        e.applied = value;
        if (!e.count) it = x.bins.erase(it); else ++it;
      }
    }
    dirty_.clear();
    marriage.resize(totals_.size());
    cohab.resize(totals_.size());
    for (size_t b = 0; b < totals_.size(); ++b) { marriage[b] = totals_[b][0]; cohab[b] = totals_[b][1]; }
  }
  // A woman's equivalent search hazard (for consistency checks).
  double hazard(int woman) const {
    if (!active(woman)) return 0;
    const Woman& x = women_[woman];
    const CellSums& c = cells_[x.cell];
    return search_->cell_rate(x.cell) * static_cast<double>(c.count) * x.weight / c.weight;
  }
  double probability(int woman) const { return probability_[woman]; }
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
  bool cohab = false, converted = false;
  double converted_from = NA_REAL;
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
                  Nullable<List> model_schedule_input = R_NilValue,
                  Nullable<List> weight_input = R_NilValue) {
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
  // Parity selects fertility cells when schedules are parity-specific; it is
  // tracked always because weight models may use it.
  std::vector<int> parity(people.size(), 0);
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
  // Weighted search under matching calibration (step 2); set with the models.
  SearchCalibration search_calibration;
  bool search_refresh_weights = false;
  // Cell-accounted fertility under birth spacing (D101).
  std::unique_ptr<FertilityCells> cells;
  // Other rate processes; search, dissolution and conversion are null when
  // not configured. One transition process per outgoing-schedule slot.
  std::unique_ptr<RateCells> mortality, search, dissolution, conversion;
  std::vector<std::unique_ptr<RateCells>> transitions;
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
        if (people[found->second].female) ++parity[found->second];
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
      if (search_refresh_weights && !fertility_driven) {
        // Weighted searches: marriage and cohabitation channels from the cells.
        std::vector<double> marriage, cohab;
        search_calibration.apply(marriage, cohab);
        if (check_cells) {
          std::vector<double> scan_m(marriage.size(), 0.0), scan_c(cohab.size(), 0.0);
          for (int w = 0; w < static_cast<int>(people.size()); ++w) {
            if (!search_calibration.active(w)) continue;
            const double h = search_calibration.hazard(w), p = search_calibration.probability(w);
            const size_t b = search_calibration.bin(w);
            scan_m[b] += h * (1 - p); scan_c[b] += h * p;
          }
          for (size_t k = 0; k < scan_m.size(); ++k)
            if (std::abs(scan_m[k] - marriage[k]) > 1e-9 * std::max(1.0, std::abs(scan_m[k])) ||
                std::abs(scan_c[k] - cohab[k]) > 1e-9 * std::max(1.0, std::abs(scan_c[k])))
              stop("Internal error: weighted search calibration differs from an explicit scan.");
        }
        calibrator->assign_channel(0, marriage);
        calibrator->assign_channel(2, cohab);
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
  // Cells stratify by group only when a process's active rates do (directly
  // or through group multipliers); otherwise groups share one set of cells.
  bool mortality_grouped = false, search_grouped = false, conversion_grouped = false;
  auto load_schedules = [&]() {
    mortality_grouped = mortality_rates.containsElementNamed(".hammel_group");
    search_grouped = union_rates.containsElementNamed(".hammel_group");
    conversion_grouped = conversion_rates.containsElementNamed(".hammel_group");
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
    if (!scheduled_models) return false;
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
    return matching_changed;
  };
  const double end = start + years;
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
  std::vector<double> tx_person, tx_time;
  std::vector<int> tx_origin, tx_dest;
  // Fertility type labels: overall, unpartnered/partnered, or three types.
  std::vector<std::string> fertility_types;
  bool fertility_grouped = false;   // fertility cells stratify by group
  // Unrealized births: (period, group, union type, parity start/end, age start/end).
  std::map<std::tuple<long long, int, std::string, double, double, double, double>, double> unrealized_counts;
  // ---- Weight models (milestone 8 step 2) ----------------------------------
  // A model gives each unit a log weight eta; selection within a chosen cell is
  // in proportion to exp(eta). Formula models are columns of factors (variable
  // code, level code or -1 for a numeric value, ancestry component); custom
  // models call an R function on the chosen cell's members.
  struct WeightSpec {
    bool custom = false, warned_na = false;
    std::string process;
    std::vector<double> coefficients;
    std::vector<std::vector<std::array<int, 3>>> factors;
    SEXP fun = R_NilValue;
  };
  std::map<std::string, WeightSpec> weight_specs;
  // Calendar schedules of weight models: intervals and their models. The
  // active model is copied into weight_specs, so references to it stay valid.
  struct WeightSchedule {
    std::vector<double> start, end;
    std::vector<WeightSpec> models;
    int active = -1;
  };
  std::map<std::string, WeightSchedule> weight_schedules;
  std::vector<std::string> weight_group_labels, weight_ancestry_labels;
  auto parse_weight_spec = [&](List x, const std::string& name) {
    WeightSpec spec;
    spec.process = name;
    spec.custom = as<bool>(x["custom"]);
    if (spec.custom) spec.fun = x["fun"];
    else {
      spec.coefficients = as<std::vector<double>>(x["coefficients"]);
      List factors = x["factors"];
      for (int k = 0; k < factors.size(); ++k) {
        IntegerMatrix m = factors[k];
        std::vector<std::array<int, 3>> column;
        for (int r = 0; r < m.nrow(); ++r) column.push_back({m(r, 0), m(r, 1), m(r, 2)});
        spec.factors.push_back(column);
      }
    }
    return spec;
  };
  auto weight_interval = [&](const WeightSchedule& w, double now) {
    for (size_t k = 0; k < w.start.size(); ++k)
      if (w.start[k] <= now && now < w.end[k]) return static_cast<int>(k);
    return -1;
  };
  if (weight_input.isNotNull()) {
    List input(weight_input.get());
    weight_group_labels = as<std::vector<std::string>>(input["group_labels"]);
    weight_ancestry_labels = as<std::vector<std::string>>(input["ancestry_labels"]);
    for (const char* name : {"mortality", "fertility", "dissolution", "conversion", "transition", "search"}) {
      if (!input.containsElementNamed(name)) continue;
      List x = input[name];
      if (x.containsElementNamed("schedule")) {
        WeightSchedule w;
        w.start = as<std::vector<double>>(x["time_start"]);
        w.end = as<std::vector<double>>(x["time_end"]);
        List models = x["models"];
        for (int k = 0; k < models.size(); ++k) w.models.push_back(parse_weight_spec(models[k], name));
        w.active = weight_interval(w, start);
        if (w.active < 0) stop(std::string("No ") + name + " weight model covers the run start.");
        // Scheduled models stay installed even when an interval is unweighted.
        weight_specs[name] = w.models[w.active];
        for (size_t k = 0; k < w.start.size(); ++k)
          for (double t : {w.start[k], w.end[k]})
            if (start < t && t < start + years) calendar_boundaries.insert(t);
        weight_schedules[name] = w;
        continue;
      }
      WeightSpec spec = parse_weight_spec(x, name);
      if (!spec.custom && spec.coefficients.empty()) continue;   // all coefficients zero: unweighted
      weight_specs[name] = spec;
    }
  }
  // Switch scheduled weight models at a calendar boundary; true if any changed.
  auto switch_weights = [&](double now) {
    bool changed = false;
    for (auto& entry : weight_schedules) {
      WeightSchedule& w = entry.second;
      const int k = weight_interval(w, now);
      if (k < 0 || k == w.active) continue;
      w.active = k;
      weight_specs[entry.first] = w.models[k];
      changed = true;
    }
    return changed;
  };
  auto weight_spec = [&](const std::string& name) -> WeightSpec* {
    auto found = weight_specs.find(name);
    return found == weight_specs.end() ? nullptr : &found->second;
  };
  auto group_label = [&](int g) -> SEXP {
    if (g < 0 || g >= static_cast<int>(weight_group_labels.size())) return NA_STRING;
    return Rf_mkChar(weight_group_labels[g].c_str());
  };
  // Variable values (codes as in R's weight_variable_codes). NA is NaN.
  // Person variables take a transition destination for group_dest.
  auto person_value = [&](int var, int component, int i, int dest) -> double {
    const Person& p = people[i];
    const int partner = current_partner[i];
    switch (var) {
      case 0: return p.female ? 1 : 0;
      case 1: return p.group;
      case 2: return partner < 0 ? 0 : (union_history[person_unions[i].back()].cohab ? 1 : 2);
      case 3: return partner >= 0 ? 1 : 0;
      case 4: return parity[i];
      case 5: return previous_union[i];
      case 6: return partner < 0 ? R_NaN : people[partner].group;
      case 7: return ancestry[static_cast<size_t>(i) * ancestry_components + component];
      case 8: return dest;
    }
    stop("Internal error: unknown person weight variable.");
  };
  auto union_value = [&](int var, int component, int u) -> double {
    const UnionRecord& x = union_history[u];
    const Person& w = people[x.wife];
    const Person& h = people[x.husband];
    switch (var) {
      case 2: return x.cohab ? 1 : 2;
      case 9: return w.group;
      case 10: return h.group;
      case 11: return parity[x.wife];
      case 12: return previous_union[x.wife];
      case 13: return previous_union[x.husband];
      case 14: return ancestry[static_cast<size_t>(x.wife) * ancestry_components + component];
      case 15: return ancestry[static_cast<size_t>(x.husband) * ancestry_components + component];
      case 16: return w.group == h.group ? 1 : 0;
      case 17: return w.born - h.born;   // husband's age minus wife's age
    }
    stop("Internal error: unknown union weight variable.");
  };
  auto warn_weight_na = [&](WeightSpec& spec, const std::string& what) {
    if (spec.warned_na) return;
    spec.warned_na = true;
    Rcpp::warning("Missing values in the " + spec.process + " weight " + what +
                  " were treated as a zero log-weight contribution.");
  };
  // A formula's log weight; a term with a missing variable contributes zero.
  auto formula_eta = [&](WeightSpec& spec, int unit, bool is_union, int dest) {
    double eta = 0;
    for (size_t k = 0; k < spec.coefficients.size(); ++k) {
      double product = 1;
      bool missing = false;
      for (const auto& f : spec.factors[k]) {
        const double v = is_union ? union_value(f[0], f[2], unit) : person_value(f[0], f[2], unit, dest);
        if (std::isnan(v)) { missing = true; break; }
        product *= f[1] == -1 ? v : (f[1] >= 0 && v == f[1] ? 1.0 : 0.0);
      }
      if (missing) { warn_weight_na(spec, "formula"); continue; }
      eta += spec.coefficients[k] * product;
    }
    return eta;
  };
  const std::array<const char*, 3> union_type_labels = {"unpartnered", "cohab", "marriage"};
  const std::array<const char*, 3> previous_labels = {"none", "separated", "widowed"};
  auto ancestry_name = [&](size_t k, const std::string& suffix) {
    const std::string& label = weight_ancestry_labels[k];
    const std::string base = label.rfind("group_", 0) == 0 ? label.substr(6) : label;
    return "ancestry_" + base + suffix;
  };
  // The members tibble a custom function receives.
  auto weight_members = [&](const std::vector<int>& units, bool is_union, int dest) {
    const int n = static_cast<int>(units.size());
    List columns;
    CharacterVector names;
    auto add = [&](const std::string& name, SEXP value) { columns.push_back(value); names.push_back(name); };
    if (!is_union) {
      NumericVector id(n); CharacterVector sex(n), group(n), type(n), previous(n), partner_group(n);
      IntegerVector partnered(n), births(n);
      for (int r = 0; r < n; ++r) {
        const int i = units[r];
        const int partner = current_partner[i];
        id[r] = people[i].id;
        sex[r] = people[i].female ? "f" : "m";
        group[r] = group_label(people[i].group);
        type[r] = union_type_labels[static_cast<int>(person_value(2, 0, i, -1))];
        partnered[r] = partner >= 0;
        births[r] = parity[i];
        previous[r] = previous_labels[previous_union[i]];
        partner_group[r] = partner < 0 ? NA_STRING : group_label(people[partner].group);
      }
      add("person_id", id); add("sex", sex); add("group", group); add("union_type", type);
      add("partnered", partnered); add("parity", births); add("previous_union", previous);
      add("partner_group", partner_group);
      for (size_t k = 0; k < ancestry_components; ++k) {
        NumericVector share(n);
        for (int r = 0; r < n; ++r) share[r] = person_value(7, static_cast<int>(k), units[r], -1);
        add(ancestry_name(k, ""), share);
      }
      if (dest >= 0) {
        CharacterVector destination(n, group_label(dest));
        add("group_dest", destination);
      }
    } else {
      NumericVector id(n), wife(n), husband(n), gap(n);
      CharacterVector group_w(n), group_h(n), previous_w(n), previous_h(n), type(n);
      IntegerVector parity_w(n), same(n);
      for (int r = 0; r < n; ++r) {
        const UnionRecord& x = union_history[units[r]];
        id[r] = x.id; wife[r] = people[x.wife].id; husband[r] = people[x.husband].id;
        group_w[r] = group_label(people[x.wife].group); group_h[r] = group_label(people[x.husband].group);
        parity_w[r] = parity[x.wife];
        previous_w[r] = previous_labels[previous_union[x.wife]];
        previous_h[r] = previous_labels[previous_union[x.husband]];
        same[r] = people[x.wife].group == people[x.husband].group;
        gap[r] = people[x.wife].born - people[x.husband].born;
        type[r] = x.cohab ? "cohab" : "marriage";
      }
      add("union_id", id); add("wife_id", wife); add("husband_id", husband);
      add("group_w", group_w); add("group_h", group_h); add("parity_w", parity_w);
      add("previous_union_w", previous_w); add("previous_union_h", previous_h);
      for (const std::string suffix : {"_w", "_h"})
        for (size_t k = 0; k < ancestry_components; ++k) {
          NumericVector share(n);
          for (int r = 0; r < n; ++r)
            share[r] = union_value(suffix == "_w" ? 14 : 15, static_cast<int>(k), units[r]);
          add(ancestry_name(k, suffix), share);
        }
      add("same_group", same); add("age_difference", gap); add("union_type", type);
    }
    columns.attr("names") = names;
    columns.attr("class") = CharacterVector::create("tbl_df", "tbl", "data.frame");
    columns.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
    return columns;
  };
  // The active rate period of a process table around now (NA when static).
  auto rate_period = [&](DataFrame rates, double now) {
    if (!rates.containsElementNamed("time_start")) return std::make_pair(NA_REAL, NA_REAL);
    NumericVector lo = rates["time_start"], hi = rates["time_end"];
    double a = R_NegInf, b = R_PosInf;
    for (int r = 0; r < lo.size(); ++r) {
      for (double t : {lo[r], hi[r]}) {
        if (t <= now) a = std::max(a, t); else b = std::min(b, t);
      }
    }
    return std::make_pair(a, b);
  };
  // Call a custom weight function and check its log weights.
  auto custom_eta = [&](WeightSpec& spec, List members, List bin, std::vector<double>& eta) {
    Function f(spec.fun);
    PutRNGstate();
    NumericVector result;
    try {
      result = as<NumericVector>(f(members, bin));
    } catch (...) {
      GetRNGstate();
      throw;
    }
    GetRNGstate();
    const int n = Rf_xlength(members[0]);
    if (result.size() != n) stop("A custom " + spec.process + " weight function returned " +
      std::to_string(result.size()) + " values for " + std::to_string(n) + " members.");
    eta.assign(n, 0.0);
    for (int r = 0; r < n; ++r) {
      if (NumericVector::is_na(result[r])) { warn_weight_na(spec, "function result"); continue; }
      if (!R_finite(result[r])) stop("A custom " + spec.process + " weight function returned an infinite value.");
      eta[r] = result[r];
    }
  };
  auto band_list = [&](List bin, bool duration, double lo, double hi, DataFrame all, double now) {
    bin[duration ? "dur_start" : "age_start"] = lo;
    bin[duration ? "dur_end" : "age_end"] = hi;
    const auto period = rate_period(all, now);
    bin["time_start"] = period.first;
    bin["time_end"] = period.second;
    return bin;
  };
  double weight_now = start;   // event time, for custom-function rate periods
  // Selection within a chosen RateCells cell for one process.
  auto make_chooser = [&](WeightSpec* spec, RateCells* process, int slot) {
    return [&, spec, process, slot](int c, const std::vector<int>& members, double uniform) -> size_t {
      const std::string& name = spec->process;
      const bool is_union = name == "dissolution" || name == "conversion";
      const int stratum = process->cell_stratum(c);
      int dest = -1;
      if (name == "transition") dest = transition_schedules[stratum][slot].dest;
      std::vector<double> eta(members.size());
      if (!spec->custom) {
        for (size_t k = 0; k < members.size(); ++k) eta[k] = formula_eta(*spec, members[k], is_union, dest);
      } else {
        List bin = List::create(_["process"] = name);
        bool duration = false;
        DataFrame all = all_mortality;
        if (name == "mortality") {
          bin["sex"] = stratum % 2 ? "f" : "m";
          bin["group"] = Rf_ScalarString(mortality_grouped ? group_label(stratum / 2) : NA_STRING);
        } else if (name == "dissolution") {
          const int pair = stratum / 2;
          bin["group_pair"] = Rf_ScalarString(dissolution_schedules.size() == 1 ? NA_STRING : Rf_mkChar(
            (weight_group_labels[pair % group_count] + ":" + weight_group_labels[pair / group_count]).c_str()));
          bin["union_type"] = stratum % 2 ? "cohab" : "marriage";
          duration = true; all = all_dissolution;
        } else if (name == "conversion") {
          bin["group"] = Rf_ScalarString(conversion_grouped ? group_label(stratum) : NA_STRING);
          duration = schedules[stratum].duration_conversion; all = all_conversion;
        } else if (name == "search") {
          const int g = previous_specific ? stratum / 3 : stratum;
          const int previous = previous_specific ? stratum % 3 : 0;
          bin["group"] = Rf_ScalarString(search_grouped ? group_label(g) : NA_STRING);
          bin["previous_union"] = Rf_ScalarString(previous_specific ? Rf_mkChar(previous_labels[previous]) : NA_STRING);
          duration = previous_specific && schedules[g].duration_search[previous]; all = all_search;
        } else {
          const auto& x = transition_schedules[stratum][slot];
          bin["group_origin"] = Rf_ScalarString(group_label(stratum));
          bin["group_dest"] = Rf_ScalarString(group_label(x.dest));
          bin["sex"] = Rf_ScalarString(x.sex < 0 ? NA_STRING : Rf_mkChar(x.sex ? "f" : "m"));
          bin["union_type"] = Rf_ScalarString(x.type.empty() ? NA_STRING : Rf_mkChar(x.type.c_str()));
          duration = transition_duration; all = all_transitions;
        }
        custom_eta(*spec, weight_members(members, is_union, dest),
          band_list(bin, duration, process->cell_start(c), process->cell_end(c), all, weight_now), eta);
      }
      return pick_log_weights(eta, uniform);
    };
  };
  // Formula models store weights (set by stored_weight at each refresh);
  // custom functions choose within the fired cell.
  auto install_weights = [&](WeightSpec* spec, RateCells* process, int slot) {
    if (spec->custom) process->set_chooser(make_chooser(spec, process, slot));
    else process->enable_stored_weights();
  };
  auto stored_weight = [&](WeightSpec* spec, RateCells* process, int unit, bool is_union, int dest) {
    if (!spec || spec->custom || !process) return;
    const double w = std::exp(formula_eta(*spec, unit, is_union, dest));
    if (!(w > 0) || !R_finite(w)) stop("A " + spec->process + " weight is not finite and positive.");
    process->set_weight(unit, w);
  };
  auto fertility_bin = [&](int c) {
    const auto info = cells->cell_info(c);
    List bin = List::create(_["process"] = "fertility");
    bin["group"] = Rf_ScalarString(fertility_grouped ? group_label(info.group) : NA_STRING);
    const std::string& type = fertility_types[info.type];
    bin["union_type"] = Rf_ScalarString(type.empty() ? NA_STRING : Rf_mkChar(type.c_str()));
    bin["parity_start"] = info.parity_start;
    bin["parity_end"] = info.parity_end;
    return band_list(bin, false, info.age_start, info.age_end, all_fertility, weight_now);
  };
  WeightSpec* fertility_weights_spec = weight_spec("fertility");
  // With matching calibration, fertility weights are stored at refresh.
  // Formula weights change only at events, so they are stored at refresh and
  // members are chosen by rejection sampling. Custom functions are evaluated when a
  // cell is chosen, except fertility under matching calibration.
  const bool fertility_refresh_weights = fertility_weights_spec &&
    (calibrator || !fertility_weights_spec->custom);
  WeightSpec* search_weights_spec = weight_spec("search");
  WeightSpec* mortality_weights_spec = weight_spec("mortality");
  WeightSpec* transition_weights_spec = weight_spec("transition");
  search_refresh_weights = search_weights_spec && calibrator && rate_initiation;
  if (search_refresh_weights && search_weights_spec->custom)
    Rcpp::warning("With weighted matching calibration, the custom search weight function is "
                  "evaluated at every refresh of every searching woman, which can be slow.");
  // A searching woman's stored weight in search cell c (calibration mode).
  std::function<double(int, int)> search_weight = [&](int woman, int c) {
    double eta;
    if (!search_weights_spec->custom) eta = formula_eta(*search_weights_spec, woman, false, -1);
    else {
      std::vector<double> values;
      std::vector<int> one = {woman};
      List bin = List::create(_["process"] = "search");
      const int stratum = search->cell_stratum(c);
      const int g = previous_specific ? stratum / 3 : stratum;
      const int previous = previous_specific ? stratum % 3 : 0;
      bin["group"] = Rf_ScalarString(search_grouped ? group_label(g) : NA_STRING);
      bin["previous_union"] = Rf_ScalarString(previous_specific ? Rf_mkChar(previous_labels[previous]) : NA_STRING);
      custom_eta(*search_weights_spec, weight_members(one, false, -1),
        band_list(bin, previous_specific && schedules[g].duration_search[previous],
                  search->cell_start(c), search->cell_end(c), all_search, weight_now), values);
      eta = values[0];
    }
    const double w = std::exp(eta);
    if (!(w > 0) || !R_finite(w)) stop("A search weight is not finite and positive.");
    return w;
  };
  if (search_refresh_weights)
    intensity_tracker->set_probability_hook([&](int woman, double p) {
      search_calibration.set_probability(woman, p);
    });
  if (fertility_refresh_weights && fertility_weights_spec->custom && calibrator)
    Rcpp::warning("With weighted matching calibration, the custom fertility weight function is "
                  "evaluated at every refresh of every woman, which can be slow.");
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
      st.stratum = cells->stratum(fertility_grouped ? p.group : 0, type, parity_specific ? parity[w] : 0);
      st.group = p.group;
      if (st.stratum < 0) stop("No fertility schedule covers a woman's current parity.");
      st.weight = heterogeneous ? fertility_weights[w] : 1;
      st.born = p.born;
      st.ready = p.last_birth + birth_interval;
      st.unpartnered = current_partner[w] < 0;
      return st;
    };
    double oldest_birth = end;
    for (const Person& person : people) oldest_birth = std::min(oldest_birth, person.born);
    fertility_grouped = fertility_rates.containsElementNamed(".hammel_group");
    const int fertility_groups = fertility_grouped ? group_count : 1;
    cells.reset(new FertilityCells(fertility_strata(fertility_rates, fertility_groups, fertility_types),
      fertility_groups, static_cast<int>(fertility_types.size()), parity_specific, end - oldest_birth,
      birth_interval, heterogeneous || fertility_refresh_weights, state, threshold, calibrator.get(),
      now, end, start));
    cells->set_pooled_groups(!fertility_grouped && group_count > 1);
    cells->set_unrealized_hook([&](int c, int group) {
      // Record one unrealized birth in its one-year period from the run start.
      long long k = static_cast<long long>(std::floor(weight_now - start));
      if (start + static_cast<double>(k + 1) <= weight_now) ++k;   // roundoff at a boundary
      const auto info = cells->cell_info(c);
      ++unrealized_counts[std::make_tuple(k, group, fertility_types[info.type], info.parity_start,
        info.parity_end, info.age_start, info.age_end)];
    });
    if (WeightSpec* spec = fertility_weights_spec) {
      if (fertility_refresh_weights) {
        cells->set_refresh_eta([&, spec](int woman, int c) {
          if (!spec->custom) return formula_eta(*spec, woman, false, -1);
          std::vector<double> eta;
          custom_eta(*spec, weight_members({woman}, false, -1), fertility_bin(c), eta);
          return eta[0];
        });
      } else {
        cells->set_select_eta([&, spec](int c, const std::vector<int>& women, std::vector<double>& eta) {
          if (spec->custom) { custom_eta(*spec, weight_members(women, false, -1), fertility_bin(c), eta); return; }
          eta.resize(women.size());
          for (size_t k = 0; k < women.size(); ++k) eta[k] = formula_eta(*spec, women[k], false, -1);
        });
      }
    }
  };
  load_cells(start);
  // Build the non-fertility rate processes from the active schedules. Strata:
  // mortality group x sex; search group (x previous-union category);
  // dissolution group pair x union type; conversion group; transition slot k
  // has one stratum per origin group (its k-th outgoing schedule).
  // The strata of each non-fertility process from the active schedules.
  struct ProcessStrata {
    std::vector<std::vector<Band>> mortality, search, dissolution, conversion;
    std::vector<std::vector<std::vector<Band>>> transitions;   // per outgoing-schedule slot
    std::vector<std::vector<std::tuple<int, int, std::string>>> transition_layout;
  };
  auto process_strata = [&]() {
    ProcessStrata x;
    for (int g = 0; g < (mortality_grouped ? group_count : 1); ++g) {
      x.mortality.push_back(schedules[g].male);
      x.mortality.push_back(schedules[g].female);
    }
    if (rate_initiation)
      for (int g = 0; g < (search_grouped ? group_count : 1); ++g) {
        if (previous_specific) for (int k = 0; k < 3; ++k) x.search.push_back(schedules[g].previous_search[k]);
        else x.search.push_back(schedules[g].search);
      }
    if (has_dissolution)
      for (size_t pair = 0; pair < dissolution_schedules.size(); ++pair) {
        x.dissolution.push_back(dissolution_schedules[pair]);
        x.dissolution.push_back(cohab_dissolution_schedules[pair]);
      }
    if (cohab_enabled)
      for (int g = 0; g < (conversion_grouped ? group_count : 1); ++g) x.conversion.push_back(schedules[g].conversion);
    size_t slots = 0;
    for (const auto& origin : transition_schedules) slots = std::max(slots, origin.size());
    for (size_t k = 0; k < slots; ++k) {
      std::vector<std::vector<Band>> strata;
      for (int g = 0; g < group_count; ++g)
        strata.push_back(k < transition_schedules[g].size() ? transition_schedules[g][k].bands :
                         std::vector<Band>());
      x.transitions.push_back(strata);
    }
    for (const auto& origin : transition_schedules) {
      std::vector<std::tuple<int, int, std::string>> layout;
      for (const auto& schedule : origin) layout.emplace_back(schedule.dest, schedule.sex, schedule.type);
      x.transition_layout.push_back(layout);
    }
    return x;
  };
  std::vector<std::vector<std::tuple<int, int, std::string>>> transition_layout;
  auto load_processes = [&](double now) {
    auto make = [&](const std::vector<std::vector<Band>>& strata) {
      return std::unique_ptr<RateCells>(new RateCells(strata, threshold, now, end));
    };
    const ProcessStrata layout = process_strata();
    transition_layout = layout.transition_layout;
    mortality = make(layout.mortality);
    search.reset();
    if (rate_initiation) {
      search = make(layout.search);
      if (search_refresh_weights) {
        search_calibration.reset(search.get(), calibrator->cells());
        // Stored weights pick the searcher; bin moves update calibration sums.
        search->enable_stored_weights();
        search->set_move_hook([&](int woman, int, int cell) {
          if (cell < 0) { search_calibration.update(woman, -1, 0, 1); return; }
          const double w = search_weights_spec->custom || !search_calibration.active(woman) ?
            search_weight(woman, cell) : search_calibration.weight(woman);
          search->set_weight(woman, w);
          search_calibration.update(woman, cell, calibrator->bin_index(people[woman].group, people[woman].born), w);
        });
      } else if (WeightSpec* spec = weight_spec("search")) {
        install_weights(spec, search.get(), -1);
      }
    }
    dissolution.reset();
    if (has_dissolution) dissolution = make(layout.dissolution);
    conversion.reset();
    if (cohab_enabled) conversion = make(layout.conversion);
    transitions.clear();
    for (size_t k = 0; k < layout.transitions.size(); ++k) {
      transitions.push_back(make(layout.transitions[k]));
      if (WeightSpec* spec = weight_spec("transition"))
        install_weights(spec, transitions.back().get(), static_cast<int>(k));
    }
    if (WeightSpec* spec = weight_spec("mortality")) install_weights(spec, mortality.get(), -1);
    if (WeightSpec* spec = weight_spec("dissolution"))
      if (dissolution) install_weights(spec, dissolution.get(), -1);
    if (WeightSpec* spec = weight_spec("conversion"))
      if (conversion) install_weights(spec, conversion.get(), -1);
  };
  // Calendar rate change: update every process's rates in place when all keep
  // their cells; otherwise return false so the caller rebuilds them.
  auto update_process_rates = [&](double now) {
    const ProcessStrata layout = process_strata();
    if (layout.transition_layout != transition_layout) return false;
    if (!mortality->same_layout(layout.mortality)) return false;
    if (search && !search->same_layout(layout.search)) return false;
    if (dissolution && !dissolution->same_layout(layout.dissolution)) return false;
    if (conversion && !conversion->same_layout(layout.conversion)) return false;
    if (layout.transitions.size() != transitions.size()) return false;
    for (size_t k = 0; k < transitions.size(); ++k)
      if (!transitions[k]->same_layout(layout.transitions[k])) return false;
    mortality->set_rates(layout.mortality, now);
    if (search) {
      search->set_rates(layout.search, now);
      if (search_refresh_weights) search_calibration.mark_all_dirty();
    }
    if (dissolution) dissolution->set_rates(layout.dissolution, now);
    if (conversion) conversion->set_rates(layout.conversion, now);
    for (size_t k = 0; k < transitions.size(); ++k) transitions[k]->set_rates(layout.transitions[k], now);
    return true;
  };
  // Recompute union u's membership in the union processes (dissolution and
  // conversion), whose units are unions identified by their history index.
  // Called when a union starts or ends and when a partner changes group.
  auto refresh_union = [&](int u, double now) {
    const UnionRecord& x = union_history[u];
    const bool ongoing = NumericVector::is_na(x.end);
    const Person& wife = people[x.wife];
    if (dissolution) {
      // Pair schedules are symmetric: both orders of a group pair share cells.
      const int low = std::min(wife.group, people[x.husband].group);
      const int high = std::max(wife.group, people[x.husband].group);
      const int pair = dissolution_schedules.size() == 1 ? 0 : low + group_count * high;
      dissolution->refresh(u, ongoing ? 2 * pair + (x.cohab ? 1 : 0) : -1, x.start, now);
      if (ongoing) stored_weight(weight_spec("dissolution"), dissolution.get(), u, true, -1);
    }
    if (conversion) {
      // Conversion schedules use the wife's group, by union duration or her age.
      const double origin = schedules[wife.group].duration_conversion ? x.start : wife.born;
      conversion->refresh(u, ongoing && x.cohab ? (conversion_grouped ? wife.group : 0) : -1, origin, now);
      if (ongoing && x.cohab) stored_weight(weight_spec("conversion"), conversion.get(), u, true, -1);
    }
  };
  // Recompute a person's membership in every rate process from their current
  // state. Called for each person whose state an event changes.
  auto refresh_person = [&](int i, double now) {
    const Person& p = people[i];
    const bool alive = NumericVector::is_na(p.died);
    const int partner = current_partner[i];
    if (p.female) cells->refresh(i, now);
    mortality->refresh(i, alive ? 2 * (mortality_grouped ? p.group : 0) + (p.female ? 1 : 0) : -1,
                       p.born, now);
    if (alive) stored_weight(mortality_weights_spec, mortality.get(), i, false, -1);
    if (search) {
      int stratum = -1;
      double origin = p.born;
      if (alive && p.female && partner < 0) {
        const int previous = previous_specific ? previous_union[i] : 0;
        const int g = search_grouped ? p.group : 0;
        stratum = previous_specific ? 3 * g + previous : g;
        if (previous_specific && schedules[p.group].duration_search[previous]) {
          origin = previous_union_end[i];
          if (!R_finite(origin)) stop("Duration search risk requires a recorded previous union ending.");
        }
      }
      search->refresh(i, stratum, origin, now);
      if (!search_refresh_weights && stratum >= 0) stored_weight(search_weights_spec, search.get(), i, false, -1);
      if (search_refresh_weights) {
        const int cell = search->cell_of(i);
        if (cell < 0) search_calibration.update(i, -1, 0, 1);
        else {
          search_calibration.set_probability(i, cohab_enabled ? current_hazard(
            &probability_schedules[p.group * 3 + previous_union[i]], p.born, now).first : 0.0);
          const double w = search_weight(i, cell);
          search->set_weight(i, w);
          search_calibration.update(i, cell, calibrator->bin_index(p.group, p.born), w);
        }
      }
    }
    if ((dissolution || conversion) && partner >= 0) refresh_union(person_unions[i].back(), now);
    if (!transitions.empty()) {
      const std::string type = partner < 0 ? "unpartnered" :
        (union_history[person_unions[i].back()].cohab ? "cohab" : "marriage");
      const auto& outgoing = transition_schedules[p.group];
      for (size_t k = 0; k < transitions.size(); ++k) {
        bool at_risk = alive && k < outgoing.size();
        if (at_risk) {
          const auto& x = outgoing[k];
          if (x.sex >= 0 && x.sex != (p.female ? 1 : 0)) at_risk = false;
          if (!x.type.empty() && x.type != type && !(x.type == "partnered" && partner >= 0)) at_risk = false;
        }
        transitions[k]->refresh(i, at_risk ? p.group : -1,
          transition_duration ? group_entry[i] : p.born, now);
        if (at_risk) stored_weight(transition_weights_spec, transitions[k].get(), i, false, outgoing[k].dest);
      }
    }
  };
  load_processes(start);
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    refresh_person(i, start);
  }
  // Event sources in same-time priority order: death, birth, dissolution,
  // conversion, search, then transition slots. Birth is the fertility cells.
  struct Source { EventType type; RateCells* cells; int slot; };
  std::vector<Source> sources;
  auto collect_sources = [&]() {
    sources.clear();
    sources.push_back({EventType::death, mortality.get(), -1});
    sources.push_back({EventType::birth, nullptr, -1});
    if (dissolution) sources.push_back({EventType::dissolution, dissolution.get(), -1});
    if (conversion) sources.push_back({EventType::conversion, conversion.get(), -1});
    if (search) sources.push_back({EventType::search, search.get(), -1});
    for (size_t k = 0; k < transitions.size(); ++k)
      sources.push_back({EventType::transition, transitions[k].get(), static_cast<int>(k)});
  };
  collect_sources();
  auto source_clock = [&](const Source& x) { return x.cells ? x.cells->next_clock() : cells->next_clock(); };
  auto source_boundary = [&](const Source& x) { return x.cells ? x.cells->next_boundary() : cells->next_boundary(); };
  double birth_count = 0, death_count = 0, dissolution_count = 0, conversion_count = 0;
  unsigned long long processed = 0;
  bool cells_consistent = true;
  // Process events in chronological order until no scheduled events remain.
  while (true) {
    // Calendar changes win exact ties with every cell/person boundary or event.
    // Conditional Exp(1) resampling is valid after survival to the boundary.
    // Candidates under the superseded rates never execute beyond it.
    if (!calendar_boundaries.empty()) {
      const double now = *calendar_boundaries.begin();
      bool first = true;
      for (const Source& x : sources)
        if (!(now <= source_boundary(x) && now <= source_clock(x))) { first = false; break; }
      if (first) {
        calendar_boundaries.erase(calendar_boundaries.begin());
        weight_now = now;
        checkUserInterrupt();
        const bool matching_changed = update_models(now);
        const bool weights_changed = switch_weights(now);
        const bool probability_changed = calendar_changes_at(all_probability, now);
        const bool rate_change = calendar_changes_at(all_mortality, now) ||
          calendar_changes_at(all_fertility, now) || calendar_changes_at(all_search, now) ||
          calendar_changes_at(all_dissolution, now) || calendar_changes_at(all_conversion, now) ||
          calendar_changes_at(all_transitions, now);
        if (!rate_change && !matching_changed && !probability_changed && !weights_changed) continue;
        for (int w = 0; w < static_cast<int>(counted_woman.size()); ++w)
          remove_counted_woman(w);
        if (probability_changed) {
          cohab_prob = calendar_slice(all_probability, now);
          load_probabilities();
        }
        if (!rate_change) {
          // Model/probability changes retain all demographic event clocks.
          // Role preferences and intensities use the newly active inputs.
          for (int w = 0; w < static_cast<int>(people.size()); ++w) {
            if (w % 1024 == 0) checkUserInterrupt();
            add_counted_woman(w, now);
          }
          // New weight models: recompute every stored weight. Intensities and
          // event clocks are unchanged, so no draws are needed.
          if (weights_changed)
            for (int i = 0; i < static_cast<int>(people.size()); ++i) {
              if (i % 1024 == 0) checkUserInterrupt();
              if (NumericVector::is_na(people[i].died)) refresh_person(i, now);
            }
          continue;
        }
        const bool fertility_change = calendar_changes_at(all_fertility, now);
        if (fertility_change) cells_consistent = cells_consistent && cells->consistent();
        for (const Source& x : sources)
          if (x.cells) cells_consistent = cells_consistent && x.cells->consistent();
        mortality_rates = calendar_slice(all_mortality, now);
        fertility_rates = calendar_slice(all_fertility, now);
        union_rates = calendar_slice(all_search, now);
        dissolution_rates = calendar_slice(all_dissolution, now);
        conversion_rates = calendar_slice(all_conversion, now);
        active_transitions = calendar_slice(all_transitions, now);
        load_schedules(); load_dissolution(); load_transitions();
        // Rates that keep every cell change in place: clocks continue and
        // nobody moves. A changed layout rebuilds that part, and its new
        // thresholds are valid after survival to the boundary (memoryless).
        bool rebuilt = false;
        if (fertility_change) {
          bool typed = false;
          if (status_fertility) {
            CharacterVector types = fertility_rates["union_type"];
            for (auto type : types) if (as<std::string>(type) == "cohab") typed = true;
          }
          const bool same = typed == typed_fertility &&
            fertility_rates.containsElementNamed(".hammel_group") == fertility_grouped &&
            cells->update_rates(fertility_strata(fertility_rates, fertility_grouped ? group_count : 1,
                                                 fertility_types), now);
          if (!same) { load_cells(now); rebuilt = true; }
        }
        if (!update_process_rates(now)) {
          load_processes(now);
          collect_sources();
          rebuilt = true;
        }
        // All table objects referenced by the tracker keep their addresses.
        // Clock selection is fixed across calendar periods by R validation.
        for (int w = 0; w < static_cast<int>(people.size()); ++w)
          add_counted_woman(w, now);
        // Rebuilt processes start empty; only the living (and their ongoing
        // unions) re-enter, so the dead need no refresh.
        if (rebuilt)
          for (int i = 0; i < static_cast<int>(people.size()); ++i) {
            if (i % 1024 == 0) checkUserInterrupt();
            if (!NumericVector::is_na(people[i].died)) continue;
            refresh_person(i, now);
          }
        continue;
      }
    }
    // Cell intensities depend on current membership, so bookkeeping boundaries
    // are processed in time order before any event at or after them. Exact
    // event ties follow the source priority order (earlier source wins).
    double boundary = R_PosInf, next_event = R_PosInf;
    size_t boundary_source = 0, event_source = 0;
    for (size_t k = 0; k < sources.size(); ++k) {
      const double b = source_boundary(sources[k]), c = source_clock(sources[k]);
      if (b < boundary) { boundary = b; boundary_source = k; }
      if (c < next_event) { next_event = c; event_source = k; }
    }
    if (boundary < end && boundary <= next_event) {
      weight_now = boundary;
      const Source& x = sources[boundary_source];
      if (x.cells) x.cells->advance_boundary(); else cells->advance_boundary();
      continue;
    }
    if (!(next_event < end)) break;
    if (++processed % 1024 == 0) checkUserInterrupt();
    const double now = next_event;
    weight_now = now;
    const Source& fired = sources[event_source];
    const EventType type = fired.type;
    // Union processes return a union index; their events act on its wife.
    const int unit = fired.cells ? fired.cells->fire_next(now) : cells->fire_next(now);
    const bool union_event = type == EventType::dissolution || type == EventType::conversion;
    if (unit < 0) continue;   // an unrealized birth: no event follows
    const int i = union_event ? union_history[unit].wife : unit;
    if (intensity_tracker) intensity_tracker->advance(now);
    if (check_cells && cells && !cells->consistent())
      stop("Internal error: fertility cell counts differ from an explicit scan.");
    if (check_cells) for (const Source& x : sources)
      if (x.cells && !x.cells->consistent())
        stop("Internal error: rate cell counts differ from an explicit scan.");
    if (type == EventType::transition) {
      const int destination = transition_schedules[people[i].group][fired.slot].dest;
      if (destination == people[i].group) stop("Invalid transition destination.");
      const int partner = current_partner[i];
      if (people[i].female) remove_counted_woman(i); else remove_male(i);
      tx_person.push_back(people[i].id); tx_time.push_back(now);
      tx_origin.push_back(people[i].group); tx_dest.push_back(destination);
      people[i].group = destination;
      group_entry[i] = now;
      if (people[i].female) add_counted_woman(i, now); else add_male(i);
      // A move changes the union's group pair and the partner's partner group.
      refresh_person(i, now);
      if (partner >= 0) refresh_person(partner, now);
    } else if (type == EventType::death) {
      // Retain the person's history, record death, and leave every process.
      people[i].died = now;
      std::vector<int> widowed, ended;
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
        widowed.push_back(i == relationship.wife ? relationship.husband : relationship.wife);
        ended.push_back(u);
      }
      if (!people[i].female) remove_male(i);
      else remove_counted_woman(i);
      refresh_person(i, now);
      for (int survivor : widowed) refresh_person(survivor, now);
      for (int u : ended) refresh_union(u, now);
      ++death_count;
    } else if (type == EventType::dissolution) {
      const int partner = current_partner[i];
      if (partner < 0) stop("Dissolution requires an ongoing union.");
      UnionRecord& relationship = union_history[unit];
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
      refresh_union(unit, now);
      refresh_person(i, now);
      refresh_person(partner, now);
    } else if (type == EventType::conversion) {
      const int partner = current_partner[i];
      const int old = unit;
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
      refresh_union(old, now);   // the cohabitation closes; the marriage enters below
      refresh_person(i, now);
      refresh_person(partner, now);
    } else if (type == EventType::search) {
      const int partner = form_union(i, now);
      refresh_person(i, now);
      if (partner >= 0) refresh_person(partner, now);
    } else {
      double father = NA_REAL;
      int partner = current_partner[i];
      if (!rate_initiation && partner < 0) {
        partner = form_union(i, now);
        if (partner >= 0) refresh_person(partner, now);
      }
      // A nonmarital father stays available: parentage alone is not a union.
      const int father_position = partner >= 0 ? partner :
        (rate_initiation ? find_candidate(i, now, father_matching, *father_draw, true, nullptr) : -1);
      if (father_position < 0 && cancel_birth_if_no_match) {
        // No child was born: retain last_birth; her union state may have changed.
        refresh_person(i, now);
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
      ++parity[i];
      counted_birth(i, now);
      // Default ffab is 0.4886 - its the signature of Berkeley training!
      bool female = R::runif(0, 1) < ffab;
      int child = static_cast<int>(people.size());
      people.push_back({++next_id, now, NA_REAL, mother, father,
                        -std::numeric_limits<double>::infinity(), female, people[i].group});
      parity.push_back(0);
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
      group_entry.push_back(now);
      ++birth_count;
      // The newborn enters its cells; the mother's spacing (and possibly her
      // union state, in on-demand mode) has changed.
      refresh_person(child, now);
      refresh_person(i, now);
    }
  }
  // Unrealized births by one-year period from the run start and cell.
  DataFrame unrealized;
  {
    const size_t n = unrealized_counts.size();
    NumericVector period_start(n), period_end(n), parity_start(n), parity_end(n),
      age_start(n), age_end(n), births(n);
    IntegerVector group(n);
    CharacterVector type(n);
    size_t r = 0;
    for (const auto& entry : unrealized_counts) {
      const auto& key = entry.first;
      const double grid_start = start + static_cast<double>(std::get<0>(key));
      period_start[r] = std::max(start, grid_start);
      period_end[r] = std::min(end, grid_start + 1);
      group[r] = std::get<1>(key);
      type[r] = std::get<2>(key);
      parity_start[r] = std::get<3>(key); parity_end[r] = std::get<4>(key);
      age_start[r] = std::get<5>(key); age_end[r] = std::get<6>(key);
      births[r] = entry.second;
      ++r;
    }
    unrealized = DataFrame::create(_["period_start"] = period_start, _["period_end"] = period_end,
      _["group"] = group, _["union_type"] = type, _["parity_start"] = parity_start,
      _["parity_end"] = parity_end, _["age_start"] = age_start, _["age_end"] = age_end,
      _["unrealized_births"] = births);
    if (cells) cells_consistent = cells_consistent && cells->consistent();
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
