// weighted_choice calibration (D066): male-group weights that make expected
// partner flows meet their targets, and the search and cohabitation
// intensities that feed it.
#ifndef HAMMEL_WEIGHT_CALIBRATION_H
#define HAMMEL_WEIGHT_CALIBRATION_H

#include <Rcpp.h>

#include <array>
#include <functional>
#include <map>
#include <vector>

#include "event_queues.h"
#include "matching.h"
#include "rate_bands.h"
#include "rate_cells.h"

namespace hammel {

using namespace Rcpp;

// Counts of unpartnered people by group and birth-time bin change only at
// events; ages are evaluated at the search time. Bin width (kWeightBinWidth)
// is a numerical tolerance, not a demographic parameter. Candidate selection
// uses exact ages.
constexpr double kZeroTargetLogWeight = -1000;  // effectively zero weight
constexpr double kMaxLogWeight = 50;            // cap for unattainable targets
// Incremental pool sums: entries below this floor are recomputed exactly at a
// search, and all sums are rebuilt after this many updates, bounding drift.
constexpr double kPoolSumFloor = 1e-6;
constexpr unsigned long long kPoolSumRefresh = 100000;

// A calibration result: log weights by matching type, and target and expected
// flow shares by weighted unit (the type for a model with ancestry endogamy,
// otherwise the category; D199). `infeasible` marks targets that no weights
// could meet (or a solve that did not converge).
struct WeightSolution {
  std::vector<double> log_weights, target, flow;
  bool infeasible = false;
};

// Unpartnered men and women by matching type and birth-time bin, women's
// search and fertility intensities by bin, and incrementally maintained pool
// sums for each registered preference table. solve() finds the weights for
// one search; channels are 0 marriage search, 1 fertility, 2 cohabitation.
// Types are matching categories unless ancestry endogamy splits them by
// profile (D199); per-type arrays grow as types are registered.
class WeightCalibrator {
  const MatchingTypes* types_;
  int groups_, categories_, bins_;
  long long first_;
  double width_, cutoff_;
  std::vector<double> men_, women_;                 // [type * bins + bin]
  std::vector<std::array<double, 3>> intensities_;  // marriage, fertility
  std::vector<std::array<int, 3>> positive_contributors_;
  // Per registered utility table: sums[(male type * C + female category) *
  // bins + female bin] of count * lag weight over included male bins, before
  // the ancestry factor. Male bins [0, included_) have reached the age
  // cutoff; older bins have lower indices.
  std::vector<const LagTable*> tables_;
  std::vector<std::vector<double>> sums_;
  int included_ = 0;
  unsigned long long updates_ = 0;
  // Men per type in included bins, and the event count at which each type
  // last went from no included men to some (held typed weights, D199).
  std::vector<double> included_men_;
  std::vector<unsigned long long> eligible_stamp_;
  unsigned long long eligibility_events_ = 0;
  void count_included(int group, double delta) {
    const bool had = included_men_[group] > 0.5;
    included_men_[group] += delta;
    if (!had && included_men_[group] > 0.5) {
      eligible_stamp_[group] = ++eligibility_events_;
    }
  }
  void spread(size_t t, int male_type, int male_bin, double delta) {
    const LagTable& table = *tables_[t];
    std::vector<double>& sums = sums_[t];
    const int male_group = types_->category(male_type);
    for (int f = 0; f < categories_; ++f) {
      // Lags for consecutive female bins are contiguous in the table.
      const double* weight =
        &table.weight[table.index(f, male_group, bins_ - 1 - male_bin)];
      double* out =
        &sums[(static_cast<size_t>(male_type) * categories_ + f) * bins_];
      for (int k = 0; k < bins_; ++k) {
        out[k] += delta * weight[k];
      }
    }
  }
  void advance(double now) {
    while (included_ < bins_ && now - centre(included_) >= cutoff_) {
      for (int g = 0; g < groups_; ++g) {
        const double count = men_[static_cast<size_t>(g) * bins_ + included_];
        if (count != 0) {
          count_included(g, count);
          for (size_t t = 0; t < tables_.size(); ++t) {
            spread(t, g, included_, count);
          }
        }
      }
      ++included_;
    }
  }
  double centre(int k) const {
    return (static_cast<double>(first_ + k) + 0.5) * width_;
  }
  size_t cell(int group, double born) const {
    const long long k =
      static_cast<long long>(std::floor(born / width_)) - first_;
    if (k < 0 || k >= bins_ || group < 0 || group >= groups_) {
      stop("Person is outside the weighted-choice calibration range.");
    }
    return static_cast<size_t>(group) * bins_ + static_cast<size_t>(k);
  }

 public:
  WeightCalibrator(
    double min_birth,
    double max_birth,
    const MatchingTypes* types,
    double width,
    double cutoff
  )
      : types_(types),
        groups_(0),
        categories_(types->categories()),
        width_(width),
        cutoff_(cutoff) {
    if (!(width > 0) || !R_finite(min_birth) || !R_finite(max_birth) ||
        max_birth < min_birth) {
      stop("Invalid weighted-choice calibration range.");
    }
    first_ = static_cast<long long>(std::floor(min_birth / width));
    const long long last =
      static_cast<long long>(std::floor(max_birth / width));
    if (last - first_ + 1 > 50000000LL / std::max(1, types->count())) {
      stop("[limits] Weighted-choice calibration range is too large.");
    }
    bins_ = static_cast<int>(last - first_ + 1);
    grow();
  }
  int bins() const { return bins_; }
  // Include the male bins that have reached the age cutoff by `now`.
  void include_through(double now) { advance(now); }
  // Typed models' held weights (D199): whether a type has men in included
  // bins, and whether it gained its first since event count `stamp`.
  unsigned long long eligibility_events() const { return eligibility_events_; }
  bool has_eligible(int type) const {
    return static_cast<size_t>(type) < included_men_.size() &&
           included_men_[type] > 0.5;
  }
  bool eligible_since(int type, unsigned long long stamp) const {
    return static_cast<size_t>(type) < eligible_stamp_.size() &&
           eligible_stamp_[type] > stamp;
  }
  // Extend per-type arrays to every registered type.
  void grow() {
    const int count = types_->count();
    if (count <= groups_) {
      return;
    }
    groups_ = count;
    const size_t n = static_cast<size_t>(count) * bins_;
    men_.resize(n, 0);
    included_men_.resize(count, 0);
    eligible_stamp_.resize(count, 0);
    women_.resize(n, 0);
    intensities_.resize(n, {0, 0, 0});
    positive_contributors_.resize(n, {0, 0, 0});
    for (auto& sums : sums_) {
      sums.resize(n * categories_, 0.0);
    }
  }
  // Register a utility table before counts are added; returns its index.
  int add_table(const LagTable* table) {
    tables_.push_back(table);
    sums_.emplace_back(static_cast<size_t>(groups_) * categories_ * bins_, 0.0);
    return static_cast<int>(tables_.size()) - 1;
  }
  // Replace a role's preference table; refresh() rebuilds sums from the
  // current pool after all roles have switched. Table storage must stay alive.
  void set_table(int id, const LagTable* table) { tables_.at(id) = table; }
  // Rebuild every pool sum from the counts, discarding accumulated roundoff.
  void refresh() {
    for (size_t t = 0; t < tables_.size(); ++t) {
      std::fill(sums_[t].begin(), sums_[t].end(), 0.0);
      for (int k = 0; k < included_; ++k) {
        for (int g = 0; g < groups_; ++g) {
          const double count = men_[static_cast<size_t>(g) * bins_ + k];
          if (count != 0) {
            spread(t, g, k, count);
          }
        }
      }
    }
    updates_ = 0;
  }
  double max_sum_difference(const WeightCalibrator& other, size_t t) const {
    double difference = 0;
    for (size_t i = 0; i < sums_[t].size(); ++i) {
      difference =
        std::max(difference, std::abs(sums_[t][i] - other.sums_[t][i]));
    }
    return difference;
  }
  void add_man(int group, double born, double delta) {
    const size_t i = cell(group, born);
    men_[i] += delta;
    const int k = static_cast<int>(i % bins_);
    if (k < included_) {
      count_included(group, delta);
      for (size_t t = 0; t < tables_.size(); ++t) {
        spread(t, group, k, delta);
      }
      if (++updates_ >= kPoolSumRefresh) {
        refresh();
      }
    }
  }
  void add_woman(int group, double born, double delta) {
    women_[cell(group, born)] += delta;
  }
  size_t bin_index(int group, double born) const { return cell(group, born); }
  size_t cells() const { return intensities_.size(); }
  // Cell-accounted fertility (D101) supplies channel totals directly before a
  // fertility-driven solve; per-woman contributions to that channel are zero.
  // Cells beyond `values` (types registered since it was filled) get zero.
  void assign_channel(int channel, const std::vector<double>& values) {
    for (size_t i = 0; i < intensities_.size(); ++i) {
      intensities_[i][channel] = i < values.size() ? values[i] : 0.0;
    }
  }
  void replace_role_intensity(
    int group,
    double born,
    const std::array<double, 3>& old,
    const std::array<double, 3>& value
  ) {
    const size_t i = cell(group, born);
    for (int c = 0; c < 3; ++c) {
      positive_contributors_[i][c] += (value[c] > 0) - (old[c] > 0);
      intensities_[i][c] += value[c] - old[c];
      // Exact zero when the last positive contribution leaves: subtraction
      // roundoff must not create phantom searchers in zero-risk cells.
      if (positive_contributors_[i][c] == 0) {
        intensities_[i][c] = 0;
      }
    }
  }
  void replace_intensity(
    int group,
    double born,
    const std::array<double, 2>& old,
    const std::array<double, 2>& value
  ) {
    replace_role_intensity(
      group, born, {old[0], old[1], 0}, {value[0], value[1], 0}
    );
  }

  // Expected male flow shares follow from each woman-bin's intensity and
  // her choice shares b[u] S[w,u] / sum_v b[v] S[w,v] over weighted units.
  // Find log weights so those shares meet the mirrored target (Newton on a
  // convex objective, warm-started; `warm` holds log weights by unit). `delta`
  // is by category. Empty or uninformative states return b = 1.
  WeightSolution solve(
    double now,
    int channel,
    int table_id,
    const std::vector<double>& delta,
    std::vector<double>& warm
  );
};

// Each counted (living, unpartnered) woman's current cohabitation probability,
// with one pending boundary per woman at her next probability band edge.
// Matching calibration splits search intensity into marriage and cohabitation
// by this probability (see SearchCalibration). Consumes no randomness.
class ProbabilityTracker {
  struct Woman {
    bool active = false;
    double born = 0;
    const std::vector<Band>* schedule = nullptr;
  };
  std::vector<Woman> women_;
  BoundaryQueue queue_;
  std::function<void(int, double)> on_change_;
  void refresh(int person, double now) {
    const Woman& w = women_[person];
    const auto p = current_hazard(w.schedule, w.born, now);
    on_change_(person, p.first);
    queue_.set(person, p.second);
  }

 public:
  explicit ProbabilityTracker(std::function<void(int, double)> on_change)
      : on_change_(std::move(on_change)) {}
  void activate(
    int person,
    double born,
    const std::vector<Band>* schedule,
    double now
  ) {
    if (person >= static_cast<int>(women_.size())) {
      women_.resize(person + 1);
    }
    women_[person] = {true, born, schedule};
    refresh(person, now);
  }
  void remove(int person) {
    if (person >= static_cast<int>(women_.size()) || !women_[person].active) {
      return;
    }
    women_[person].active = false;
    queue_.clear(person);
  }
  void advance(double now) {
    while (queue_.top() <= now) {
      const auto entry = queue_.pop();
      refresh(entry.second, entry.first);
    }
  }
};

// Weighted search intensities for matching calibration (D129). With search
// weights, a woman's search hazard is rate_c * N_c * w_i / W_c, so
// her calibration contribution depends on her whole cell. Each search cell
// keeps, per calibration bin, the member count and the sums of w (1 - p) and
// w p (p: her cohabitation probability). A bin's marriage and cohabitation
// intensities are sum over cells of rate_c * N_c / W_c times those sums.
// Cells whose contents changed are recomputed before the next search solve.
class SearchCalibration {
  struct Woman {
    int cell = -1;
    size_t bin = 0;
    double weight = 1;
    bool active = false;
  };
  // Per calibration bin within a cell: members, sums, and the last applied
  // contribution (kept until applied as zero, so stale entries are removed).
  struct BinSums {
    int count = 0;
    double marriage = 0, cohab = 0;
    std::array<double, 2> applied = {0, 0};
  };
  struct CellSums {
    int count = 0;
    double weight = 0;
    std::map<size_t, BinSums> bins;
    bool dirty = false;
  };
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
    if (c >= static_cast<int>(cells_.size())) {
      cells_.resize(c + 1);
    }
    return cells_[c];
  }
  void touch(int c) {
    CellSums& x = cell(c);
    if (!x.dirty) {
      x.dirty = true;
      dirty_.push_back(c);
    }
  }
  void remove(int woman) {
    Woman& x = women_[woman];
    if (!x.active) {
      return;
    }
    CellSums& c = cell(x.cell);
    const double p = probability_[woman];
    BinSums& b = c.bins[x.bin];
    b.marriage -= x.weight * (1 - p);
    b.cohab -= x.weight * p;
    if (--b.count == 0) {
      b.marriage = b.cohab = 0;
    }
    c.weight -= x.weight;
    if (--c.count == 0) {
      c.weight = 0;
    }
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
    for (Woman& x : women_) {
      x.active = false;
    }
    cells_.clear();
    dirty_.clear();
    totals_.assign(bins, {0, 0});
    positive_.assign(bins, {0, 0});
  }
  // Place a woman in a search cell (-1: not searching) with her bin and weight.
  void update(int woman, int c, size_t bin, double weight) {
    grow(woman);
    if (c >= 0 && bin >= totals_.size()) {
      totals_.resize(bin + 1, {0, 0});
      positive_.resize(bin + 1, {0, 0});
    }
    const Woman& x0 = women_[woman];
    if (x0.active ? (c == x0.cell && bin == x0.bin && weight == x0.weight)
                  : c < 0) {
      return;
    }
    remove(woman);
    if (c < 0) {
      return;
    }
    Woman& x = women_[woman];
    x.cell = c;
    x.bin = bin;
    x.weight = weight;
    add(woman);
  }
  void set_probability(int woman, double p) {
    grow(woman);
    if (probability_[woman] == p) {
      return;
    }
    const bool active = women_[woman].active;
    remove(woman);
    probability_[woman] = p;
    if (active) {
      add(woman);
    }
  }
  double weight(int woman) const { return women_[woman].weight; }
  // Search rates changed in place: every cell's contribution is recomputed.
  void mark_all_dirty() {
    for (size_t c = 0; c < cells_.size(); ++c) {
      touch(static_cast<int>(c));
    }
  }
  size_t bin(int woman) const { return women_[woman].bin; }
  bool active(int woman) const {
    return woman < static_cast<int>(women_.size()) && women_[woman].active;
  }
  // Recompute changed cells and return marriage and cohabitation totals by bin.
  void apply(std::vector<double>& marriage, std::vector<double>& cohab) {
    for (int c : dirty_) {
      CellSums& x = cells_[c];
      x.dirty = false;
      const double scale =
        x.count > 0 && x.weight > 0
          ? search_->cell_rate(c) * static_cast<double>(x.count) / x.weight
          : 0.0;
      for (auto it = x.bins.begin(); it != x.bins.end();) {
        BinSums& e = it->second;
        const size_t b = it->first;
        const std::array<double, 2> value = {
          e.count ? scale * std::max(0.0, e.marriage) : 0.0,
          e.count ? scale * std::max(0.0, e.cohab) : 0.0};
        for (int k = 0; k < 2; ++k) {
          positive_[b][k] += (value[k] > 0) - (e.applied[k] > 0);
          totals_[b][k] += value[k] - e.applied[k];
          // Exact zero when no cell contributes: no phantom roundoff intensity.
          if (positive_[b][k] == 0) {
            totals_[b][k] = 0;
          }
        }
        e.applied = value;
        if (!e.count) {
          it = x.bins.erase(it);
        } else {
          ++it;
        }
      }
    }
    dirty_.clear();
    marriage.resize(totals_.size());
    cohab.resize(totals_.size());
    for (size_t b = 0; b < totals_.size(); ++b) {
      marriage[b] = totals_[b][0];
      cohab[b] = totals_[b][1];
    }
  }
  // A woman's equivalent search hazard (for consistency checks).
  double hazard(int woman) const {
    if (!active(woman)) {
      return 0;
    }
    const Woman& x = women_[woman];
    const CellSums& c = cells_[x.cell];
    return search_->cell_rate(x.cell) * static_cast<double>(c.count) *
           x.weight / c.weight;
  }
};

}  // namespace hammel

#endif  // HAMMEL_WEIGHT_CALIBRATION_H
