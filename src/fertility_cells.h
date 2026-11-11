// Cell-accounted fertility under birth spacing (D101) with persistent
// fertility heterogeneity (hetfert) and weight models (D129).
#ifndef HAMMEL_FERTILITY_CELLS_H
#define HAMMEL_FERTILITY_CELLS_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "event_queues.h"
#include "rate_bands.h"
#include "weight_calibration.h"

namespace hammel {

using namespace Rcpp;

// A cell is one age sub-band (at most one year, split at whole ages and input
// band edges) of one effective fertility schedule (group, union type, and
// actual parity when schedules are parity-specific; pooled bands are split by
// actual parity). Each cell generates births at m_c * N_c, counting
// all living members including those inside the birth interval, and assigns
// each birth to an eligible member, weighted by persistent hetfert when
// enabled (selection probability, not a rate multiplier). Without hetfert,
// selection is uniform and no weight tree is allocated. A cell with members but
// no eligible member keeps its intensity: when it is drawn, no child is born
// and one unrealized birth is recorded. Clocks keep an Exp(1) threshold
// consumed at the current intensity, so count changes recompute firing times
// without drawing random numbers.
//
// Upper-open fertility bands stay one cell when birth_interval * rate is below
// this share (a numerical tolerance, not a demographic parameter).
constexpr double kOpenBandSplitTolerance = 0.001;

// Log-space form of the accepted latent-normal inheritance; valid also at
// sigma=0 for supplied nonunit weights. Endpoints avoid needless RNG draws.
double inherit_fertility_weight(double maternal, double variance, double alpha);

// Compensated bin sums retain small weights when a dominant woman leaves.
struct FertilityWeightSum {
  double sum = 0, correction = 0;
  void add(double value) {
    const double next = sum + value;
    correction += std::abs(sum) >= std::abs(value) ? (sum - next) + value
                                                   : (value - next) + sum;
    sum = next;
  }
  double value() const { return sum + correction; }
};

// The fertility cells of one run, with the single birth clock. Engine state
// reaches it through the `state` callback at each refresh.
class FertilityCells {
 public:
  // A woman's engine state as fertility sees it.
  struct State {
    int stratum = -1;  // -1: not a living woman
    int group = 0;     // her matching group (calibration bins)
    int report = 0;    // her categories in all group dimensions, combined
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
    std::vector<int> eligible;  // swap-remove list of women
    FertilityWeightTree weights;
    double lambda = 0;
    std::map<size_t, std::pair<int, FertilityWeightSum>>
      counted;  // count and weight per unpartnered calibration bin
    std::map<int, int> by_group;  // members by reporting categories
                                  // (attributes unrealized births)
  };
  struct Member {
    int cell = -1, slot = -1, report = 0;
    bool eligible = false, counted = false;
    size_t bin = 0;
    double weight = 1;
  };
  std::vector<Source> sources_;
  std::vector<Stratum> strata_;
  std::vector<Cell> cells_;
  std::vector<Member> members_;
  std::vector<State> states_;  // engine state at each woman's last refresh
  // Weight models (D129). Without matching calibration, log weights are
  // computed for the chosen cell's eligible women at selection. With
  // calibration, refresh_eta_ folds them into each woman's stored weight.
  std::function<void(int, const std::vector<int>&, std::vector<double>&)>
    select_eta_;
  std::function<double(int, int)> refresh_eta_;
  // lookup_[group * types + type][actual parity] -> stratum (-1 until created);
  // without parity-specific schedules, index 0 serves every parity.
  int types_;
  bool split_parity_, weighted_;
  bool pooled_groups_ = false;  // rates shared by groups: one set of cells
  double max_age_, interval_;
  std::vector<std::vector<int>> lookup_;
  BoundaryQueue boundaries_;
  // One birth clock for all cells; every cell contributes rate * members.
  IntensityTree intensity_;
  ProcessClock clock_;
  std::function<State(int)> state_;
  std::function<double()> draw_threshold_;
  std::vector<double> eta_;  // scratch log weights for selection
  const WeightCalibrator* calibrator_;
  double end_;
  // Called with (cell, group) when a cell without eligible women is drawn.
  std::function<void(int, int)> on_unrealized_;
  unsigned long long updates_ = 0;

  // Recompute the cell's intensity from current counts, first spending the
  // clock for the interval since the last change.
  void touch(int c, double now) {
    Cell& x = cells_[c];
    x.lambda = x.rate * x.members;
    clock_.spend(intensity_, now);
    intensity_.set(static_cast<size_t>(c), x.lambda);
    clock_.arm(intensity_, draw_threshold_);
  }
  void leave(int woman, int& changed) {
    Member& m = members_[woman];
    if (m.cell < 0) {
      return;
    }
    Cell& x = cells_[m.cell];
    --x.members;
    if (pooled_groups_ && --x.by_group[m.report] == 0) {
      x.by_group.erase(m.report);
    }
    if (m.eligible) {
      const int last = x.eligible.back();
      x.eligible[m.slot] = last;
      members_[last].slot = m.slot;
      x.eligible.pop_back();
      if (weighted_) {
        x.weights.remove(m.slot);
      }
    }
    if (m.counted) {
      auto& bin = x.counted[m.bin];
      if (weighted_) {
        bin.second.add(-m.weight);
      }
      if (--bin.first == 0) {
        x.counted.erase(m.bin);
      }
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
  static std::vector<Band>
  split_ages(const std::vector<Band>& bands, double max_age, double interval) {
    std::vector<Band> out;
    for (const Band& band : bands) {
      // Zero-rate bands generate no births and need no members: treat as gaps.
      if (!(band.rate > 0)) {
        continue;
      }
      double lo = band.start;
      if (!std::isfinite(band.end) &&
          interval * band.rate < kOpenBandSplitTolerance) {
        out.push_back(band);
        continue;
      }
      const double cap =
        std::isfinite(band.end) ? band.end : std::max(lo, std::ceil(max_age));
      while (lo < cap) {
        const double next = std::min(cap, std::floor(lo) + 1);
        out.push_back({lo, next, band.rate});
        lo = next;
      }
      if (!std::isfinite(band.end)) {
        out.push_back({lo, band.end, band.rate});
      }
    }
    return out;
  }
  int add_stratum(
    const Source& source,
    double parity_start,
    double parity_end
  ) {
    Stratum st;
    st.group = source.group;
    st.type = source.type;
    st.parity_start = parity_start;
    st.parity_end = parity_end;
    st.bands = split_ages(source.bands, max_age_, interval_);
    st.first_cell = static_cast<int>(cells_.size());
    const int k = static_cast<int>(strata_.size());
    for (const Band& band : st.bands) {
      Cell x;
      x.stratum = k;
      x.rate = band.rate;
      x.age_start = band.start;
      x.age_end = band.end;
      cells_.push_back(x);
    }
    strata_.push_back(st);
    return k;
  }

  FertilityCells(
    std::vector<Source> sources,
    int groups,
    int types,
    bool split_parity,
    double max_age,
    double interval,
    bool weighted,
    std::function<State(int)> state,
    std::function<double()> draw_threshold,
    const WeightCalibrator* calibrator,
    double start,
    double end
  )
      : sources_(std::move(sources)),
        types_(types),
        split_parity_(split_parity),
        weighted_(weighted),
        max_age_(max_age),
        interval_(interval),
        lookup_(static_cast<size_t>(groups) * types),
        state_(std::move(state)),
        draw_threshold_(std::move(draw_threshold)),
        calibrator_(calibrator),
        end_(end) {
    clock_.start(start);
    if (!split_parity_) {
      for (const Source& source : sources_) {
        lookup_[static_cast<size_t>(source.group) * types_ + source.type] = {
          add_stratum(source, source.parity_start, source.parity_end)};
      }
    }
  }
  // Input parity bands partition [0, Inf) within each group/type (validated in
  // R). Actual-parity strata are created when a woman first reaches a parity.
  int stratum(int group, int type, int parity) {
    auto& list = lookup_[static_cast<size_t>(group) * types_ + type];
    if (!split_parity_) {
      return list.empty() ? -1 : list[0];
    }
    if (parity < 0) {
      return -1;
    }
    if (parity < static_cast<int>(list.size()) && list[parity] >= 0) {
      return list[parity];
    }
    for (const Source& source : sources_) {
      if (source.group != group || source.type != type ||
          !(source.parity_start <= parity && parity < source.parity_end)) {
        continue;
      }
      if (parity >= static_cast<int>(list.size())) {
        list.resize(parity + 1, -1);
      }
      list[parity] = add_stratum(source, parity, parity + 1.0);
      return list[parity];
    }
    return -1;
  }
  // Recompute a woman's cell, eligibility and calibration bin from engine
  // state.
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
      const Stratum& st = strata_[s.stratum];
      const int band = locate_band(st.bands, s.born, now, next);
      if (band >= 0) {
        cell = st.first_cell + band;
      }
    }
    const bool eligible = cell >= 0 && now >= s.ready;
    if (cell >= 0 && !eligible) {
      next = std::min(next, s.ready);
    }
    const bool counted = calibrator_ && eligible && s.unpartnered;
    const size_t bin = counted ? calibrator_->bin_index(s.group, s.born) : 0;
    double weight = s.weight;
    if (refresh_eta_ && cell >= 0) {
      weight *= std::exp(refresh_eta_(woman, cell));
      if (!(weight > 0) || !R_finite(weight)) {
        stop("A fertility weight is not finite and positive.");
      }
    }
    Member& m = members_[woman];
    if (cell >= 0 && m.cell == cell && m.report != s.report && pooled_groups_) {
      // Only her reporting categories changed (a group transition): update the
      // cell's member counts in place, keeping her place among eligible women.
      Cell& x = cells_[cell];
      if (--x.by_group[m.report] == 0) {
        x.by_group.erase(m.report);
      }
      ++x.by_group[s.report];
      m.report = s.report;
    }
    if (m.cell != cell || m.eligible != eligible || m.counted != counted ||
        (counted && m.bin != bin) || (cell >= 0 && weight != m.weight)) {
      int old_cell = -1;
      leave(woman, old_cell);
      if (cell >= 0) {
        Cell& x = cells_[cell];
        ++x.members;
        if (pooled_groups_) {
          ++x.by_group[s.report];
        }
        m.cell = cell;
        m.eligible = eligible;
        m.counted = counted;
        m.bin = bin;
        m.weight = weight;
        m.report = s.report;
        if (eligible) {
          m.slot = static_cast<int>(x.eligible.size());
          x.eligible.push_back(woman);
          if (weighted_) {
            x.weights.add(m.weight);
          }
        }
        if (counted) {
          ++x.counted[bin].first;
          if (weighted_) {
            x.counted[bin].second.add(m.weight);
          }
        }
      }
      if (old_cell >= 0) {
        touch(old_cell, now);
      }
      if (cell >= 0 && cell != old_cell) {
        touch(cell, now);
      }
    }
    boundaries_.set(woman, next);
  }
  double next_boundary() { return boundaries_.top(); }
  void advance_boundary() {
    const auto entry = boundaries_.pop();
    place(entry.second, entry.first);
    if (++updates_ % 1024 == 0) {
      checkUserInterrupt();
    }
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
          if ((target -= g.second) < 0) {
            break;
          }
        }
      }
      clock_.fired(now);
      clock_.arm(intensity_, draw_threshold_);
      if (on_unrealized_) {
        on_unrealized_(c, group);
      }
      return -1;
    }
    size_t pick = weighted_ ? x.weights.choose(uniform)
                            : uniform_index(uniform, x.eligible.size());
    if (select_eta_) {
      // Equal formula weights keep the pick above (hetfert or uniform).
      select_eta_(c, x.eligible, eta_);
      bool equal = true;
      for (double e : eta_) {
        if (e != eta_[0]) {
          equal = false;
          break;
        }
      }
      if (!equal) {
        if (weighted_) {
          for (size_t k = 0; k < eta_.size(); ++k) {
            eta_[k] += std::log(members_[x.eligible[k]].weight);
          }
        }
        pick = pick_log_weights(eta_, uniform);
      }
    }
    const int woman = x.eligible[pick];
    clock_.fired(now);
    clock_.arm(intensity_, draw_threshold_);
    return woman;
  }
  void set_select_eta(
    std::function<void(int, const std::vector<int>&, std::vector<double>&)> f
  ) {
    select_eta_ = std::move(f);
  }
  void set_refresh_eta(std::function<double(int, int)> f) {
    refresh_eta_ = std::move(f);
  }
  void set_pooled_groups(bool pooled) { pooled_groups_ = pooled; }
  // A cell's group, union-type index, parity range and age band.
  struct CellInfo {
    int group, type;
    double parity_start, parity_end, age_start, age_end;
  };
  CellInfo cell_info(int c) const {
    const Cell& x = cells_[c];
    const Stratum& st = strata_[x.stratum];
    return {st.group,      st.type,     st.parity_start,
            st.parity_end, x.age_start, x.age_end};
  }
  // Calibration fertility totals: each unpartnered eligible member contributes
  // her equivalent weighted hazard in her calibration bin.
  std::vector<double> fill_calibration(WeightCalibrator& calibrator) const {
    std::vector<double> values(calibrator.cells(), 0.0);
    for (const Cell& x : cells_) {
      if (!(x.lambda > 0) || x.eligible.empty() || x.counted.empty()) {
        continue;
      }
      if (weighted_) {
        const double total = x.weights.total();
        for (const auto& bin : x.counted) {
          values[bin.first] += x.lambda * (bin.second.second.value() / total);
        }
      } else {
        const double hazard = x.lambda / static_cast<double>(x.eligible.size());
        for (const auto& bin : x.counted) {
          values[bin.first] += hazard * bin.second.first;
        }
      }
    }
    calibrator.assign_channel(1, values);
    return values;
  }
  // A woman's equivalent individual hazard (D101 rule 4), zero inside
  // spacing.
  double hazard(int woman) const {
    if (woman >= static_cast<int>(members_.size())) {
      return 0;
    }
    const Member& m = members_[woman];
    if (m.cell < 0 || !m.eligible) {
      return 0;
    }
    const Cell& x = cells_[m.cell];
    return weighted_ ? x.lambda * (m.weight / x.weights.total())
                     : x.lambda / static_cast<double>(x.eligible.size());
  }
  void set_unrealized_hook(std::function<void(int, int)> f) {
    on_unrealized_ = std::move(f);
  }
  // Calendar rate changes: update rates in place when the new sources keep
  // every stratum's cells (same sources and, after age splitting, the same
  // band edges). Returns false, changing nothing, otherwise.
  bool update_rates(const std::vector<Source>& sources, double now) {
    if (sources.size() != sources_.size()) {
      return false;
    }
    for (size_t i = 0; i < sources.size(); ++i) {
      if (sources[i].group != sources_[i].group ||
          sources[i].type != sources_[i].type ||
          sources[i].parity_start != sources_[i].parity_start ||
          sources[i].parity_end != sources_[i].parity_end) {
        return false;
      }
    }
    std::vector<std::vector<Band>> staged(strata_.size());
    for (size_t k = 0; k < strata_.size(); ++k) {
      const Stratum& st = strata_[k];
      const Source* source = nullptr;
      for (const Source& x : sources) {
        if (x.group == st.group && x.type == st.type &&
            x.parity_start <= st.parity_start &&
            st.parity_start < x.parity_end) {
          source = &x;
          break;
        }
      }
      if (!source) {
        return false;
      }
      staged[k] = split_ages(source->bands, max_age_, interval_);
      if (staged[k].size() != st.bands.size()) {
        return false;
      }
      for (size_t i = 0; i < staged[k].size(); ++i) {
        if (staged[k][i].start != st.bands[i].start ||
            staged[k][i].end != st.bands[i].end) {
          return false;
        }
      }
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
      if (m.cell < 0) {
        continue;
      }
      ++members[m.cell];
      if (m.eligible) {
        ++eligible[m.cell];
        weights[m.cell] += m.weight;
        if (cells_[m.cell].eligible[m.slot] != static_cast<int>(w)) {
          return false;
        }
      }
    }
    for (size_t c = 0; c < cells_.size(); ++c) {
      if (members[c] != cells_[c].members ||
          eligible[c] != static_cast<int>(cells_[c].eligible.size())) {
        return false;
      }
      if (intensity_.value(c) != cells_[c].lambda) {
        return false;
      }
    }
    if (weighted_) {
      for (size_t c = 0; c < cells_.size(); ++c) {
        if (std::abs(weights[c] - cells_[c].weights.total()) >
            1e-9 * std::max(weights[c], 1.0)) {
          return false;
        }
      }
    }
    return true;
  }
};

// Input schedules for cells from the effective (group-multiplied) fertility
// table.
std::vector<FertilityCells::Source> fertility_strata(
  DataFrame rates,
  int groups,
  const std::vector<std::string>& type_labels
);

}  // namespace hammel

#endif  // HAMMEL_FERTILITY_CELLS_H
