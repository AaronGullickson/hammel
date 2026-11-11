// Partner matching: candidate utilities, the lag-grid utility table, the
// pool of unpartnered men, and exact full-pool conditional-logit draws.
#ifndef HAMMEL_MATCHING_H
#define HAMMEL_MATCHING_H

#include <Rcpp.h>

#include <functional>
#include <vector>

namespace hammel {

using namespace Rcpp;

// Width, in years, of the birth-time bins that group men for exact partner
// draws and weighted-choice calibration (D066). A numerical tolerance, not
// a demographic parameter: candidate selection always uses exact ages.
constexpr double kWeightBinWidth = 0.25;

// Evaluate candidate utilities without rounding age differences. Common linear
// and raw quadratic formulas use an exact native path; other formulas call the
// R model-matrix evaluator with candidate data only, never mutable engine
// state.
NumericVector matching_utilities(
  const NumericVector& gaps,
  const List& matching,
  int female_group = -1,
  IntegerVector male_groups = IntegerVector(),
  IntegerVector female_groups = IntegerVector()
);

// Pair utilities on the fixed grid of bin-centre lags. Lag (female bin - male
// bin) * width is male age minus female age, so the grid never changes during
// a run and custom R formulas are evaluated once. Excluded pairs are -Inf.
struct LagTable {
  int groups = 0, span = 0;
  std::vector<double> utility,
    weight;  // weight = exp(utility - global maximum)
  // allowed[female * groups + male]: the pair is not prohibited. block[g]:
  // the separate market (D138) of group g; groups connected by allowed pairs
  // share a block.
  std::vector<char> allowed;
  std::vector<int> block;
  int blocks = 1;
  size_t index(int female, int male, int lag) const {
    return (static_cast<size_t>(female) * groups + male) * span + lag;
  }
};
LagTable build_lag_table(const List& model, int groups, int bins, double width);

// Exact full-pool partner draws (D072). Living, unpartnered men are kept in
// cells by group and birth-time bin. A draw picks a cell with probability
// proportional to its size times the largest weight any member could have for
// this woman, picks a member uniformly, accepts him with probability exp(U -
// cell maximum), and screens only him. This samples conditional-logit choice
// over the full eligible pool exactly, without a candidate sample.
class MalePool {
  int groups_ = 0, bins_ = 0;
  long long first_ = 0;
  double width_ = 0;
  std::vector<std::vector<int>> cells_;
  std::vector<int> slot_, cell_;
  size_t size_ = 0;

 public:
  MalePool() = default;
  MalePool(double min_birth, double max_birth, int groups, double width)
      : groups_(groups), width_(width) {
    if (!(width > 0) || !R_finite(min_birth) || !R_finite(max_birth) ||
        max_birth < min_birth || groups < 1) {
      stop("Invalid male pool range.");
    }
    first_ = static_cast<long long>(std::floor(min_birth / width));
    const long long last =
      static_cast<long long>(std::floor(max_birth / width));
    if (last - first_ + 1 > 50000000LL / groups) {
      stop("Male pool birth-time range is too large.");
    }
    bins_ = static_cast<int>(last - first_ + 1);
    cells_.assign(static_cast<size_t>(groups) * bins_, std::vector<int>());
  }
  int groups() const { return groups_; }
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
    if (slot_[person] >= 0) {
      return;
    }
    const long long k =
      static_cast<long long>(std::floor(born / width_)) - first_;
    if (k < 0 || k >= bins_ || group < 0 || group >= groups_) {
      stop("Person is outside the male pool range.");
    }
    const int c = group * bins_ + static_cast<int>(k);
    slot_[person] = static_cast<int>(cells_[c].size());
    cell_[person] = c;
    cells_[c].push_back(person);
    ++size_;
  }
  void remove(int person) {
    if (!contains(person)) {
      return;
    }
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
    const double k =
      std::floor((now - cutoff) / width_) - static_cast<double>(first_);
    return static_cast<int>(
      std::min(static_cast<double>(bins_ - 1), std::max(-1.0, k))
    );
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
  std::vector<double> grid;  // [(female * groups + male) * grid_size + j]

  double offset(int f, int m) const {
    return offsets.nrow() ? offsets(f, m) : 0.0;
  }
  bool allowed(int f, int m) const {
    return !(prevent_all && f != m) && !(blocked.nrow() && blocked(f, m));
  }
  double native_value(double gap) const {
    return linear * gap + quadratic * gap * gap;
  }
  // Exact utilities, including pair effects, for one female/male group pair.
  NumericVector utilities(int f, int m, const NumericVector& gaps) const {
    if (native) {
      NumericVector out(gaps.size());
      for (R_xlen_t j = 0; j < gaps.size(); ++j) {
        out[j] = native_value(gaps[j]) + offset(f, m);
        if (!R_finite(out[j])) {
          stop("Matching ages and utilities must be finite.");
        }
      }
      return out;
    }
    return matching_utilities(
      gaps, model, -1, IntegerVector(gaps.size(), m),
      IntegerVector(gaps.size(), f)
    );
  }
  // The exact utility of one candidate; native models avoid allocating.
  double utility(int f, int m, double gap) const {
    if (!native) {
      return utilities(f, m, NumericVector::create(gap))[0];
    }
    const double value = native_value(gap) + offset(f, m);
    if (!R_finite(value)) {
      stop("Matching ages and utilities must be finite.");
    }
    return value;
  }
  // Upper bound on utility for male-minus-female age gaps in [lo, hi].
  double bound(int f, int m, double lo, double hi) const {
    if (native) {
      double best = std::max(native_value(lo), native_value(hi));
      if (quadratic < 0) {
        const double vertex = -linear / (2 * quadratic);
        if (vertex > lo && vertex < hi) {
          best = native_value(vertex);
        }
      }
      return best + offset(f, m);
    }
    const int j0 =
      std::max(0, static_cast<int>(std::floor((lo - grid_start) / grid_step)));
    const int j1 = std::min(
      grid_size - 1, static_cast<int>(std::ceil((hi - grid_start) / grid_step))
    );
    const double* values =
      &grid[(static_cast<size_t>(f) * groups + m) * grid_size];
    double best = R_NegInf;
    for (int j = j0; j <= j1; ++j) {
      best = std::max(best, values[j]);
    }
    return best;
  }
};
DrawModel make_draw_model(const List& model, int groups, double gap_range);

// Reusable scratch storage for draw_partner(), one per run.
struct DrawWorkspace {
  std::vector<int> group, bin;
  std::vector<double> top, cumulative;
};

// Draw one partner for a woman of group `female_group` born `female_born`.
// Returns a person index, or -1 when no eligible man exists. `log_weights`
// holds one log group weight per group (zeros for female choice).
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
);

}  // namespace hammel

#endif  // HAMMEL_MATCHING_H
