// Cell-scheduled rate processes (D128): mortality, union search,
// dissolution, conversion and group transitions.
#ifndef HAMMEL_RATE_CELLS_H
#define HAMMEL_RATE_CELLS_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "event_queues.h"
#include "rate_bands.h"

namespace hammel {

using namespace Rcpp;

// One demographic rate process (mortality, union search, dissolution,
// conversion, or one group-transition slot). Each at-risk unit (a person, or a
// union identified by its history index) belongs to one cell: one nonzero
// input band of its schedule (stratum), measured from the unit's own clock
// origin (birth, union start, last union end, or category entry). With no
// eligibility rules every member's hazard is the cell rate, so splitting a
// band could not change results, only add boundary updates. Zero-rate bands
// are gaps. A cell generates events at rate * members; one process clock fires
// at the total intensity, the cell is chosen in proportion to its intensity,
// and the unit uniformly within it, or by its weight under a weight model
// (D129).
class RateCells {
  struct Cell {
    int stratum = 0;
    double rate = 0, start = 0, end = 0;  // band edges on the unit's clock
    std::vector<int> members;             // swap-remove list of units
    double weight_bound = 0;  // upper bound on members' stored weights
  };
  struct Member {
    int cell = -1, slot = -1, stratum = -1;
    double origin = 0, weight = 1;
  };
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
  // Custom weight-function selection: returns the index of the chosen member
  // of a cell's member list, given the cell and one uniform draw.
  std::function<size_t(int, const std::vector<int>&, double)> chooser_;
  // Called after a unit changes cell (old, new; -1 for none).
  std::function<void(int, int, int)> on_move_;
  double end_;
  unsigned long long updates_ = 0;
  void update(int c, double now) {
    clock_.spend(intensity_, now);
    intensity_.set(
      static_cast<size_t>(c),
      cells_[c].rate * static_cast<double>(cells_[c].members.size())
    );
    clock_.arm(intensity_, draw_threshold_);
  }

 public:
  RateCells(
    const std::vector<std::vector<Band>>& schedules,
    std::function<double()> draw_threshold,
    double start,
    double end
  )
      : draw_threshold_(std::move(draw_threshold)), end_(end) {
    for (const auto& bands : schedules) {
      std::vector<Band> cells;
      for (const Band& band : bands) {
        if (band.rate > 0) {
          cells.push_back(band);
        }
      }
      strata_.push_back(std::move(cells));
      first_.push_back(static_cast<int>(cells_.size()));
      for (const Band& band : strata_.back()) {
        Cell x;
        x.stratum = static_cast<int>(strata_.size()) - 1;
        x.rate = band.rate;
        x.start = band.start;
        x.end = band.end;
        cells_.push_back(std::move(x));
      }
    }
    if (!cells_.empty()) {
      intensity_.set(cells_.size() - 1, 0);
    }
    clock_.start(start);
  }
  // Place a unit in the cell of schedule `stratum` (-1: not at risk) that
  // holds its clock value now - origin, and queue its next band boundary.
  void refresh(int person, int stratum, double origin, double now) {
    if (person >= static_cast<int>(members_.size())) {
      members_.resize(person + 1);
    }
    if (stratum >= 0 && !R_finite(origin)) {
      stop("A rate process requires a finite clock origin.");
    }
    // Fast path: most refreshes follow events that leave this process's state
    // unchanged. Same stratum and origin, still inside the current cell's band:
    // the cell and the queued boundary are already correct.
    Member& m = members_[person];
    if (stratum == m.stratum &&
        (stratum < 0 ? m.cell < 0
                     : origin == m.origin && m.cell >= 0 &&
                         origin + cells_[m.cell].start <= now &&
                         now < origin + cells_[m.cell].end)) {
      return;
    }
    int cell = -1;
    double next = R_PosInf;
    if (stratum >= 0) {
      const int band = locate_band(strata_[stratum], origin, now, next);
      if (band >= 0) {
        cell = first_[stratum] + band;
      }
    }
    m.stratum = stratum;
    m.origin = origin;
    if (m.cell != cell) {
      const int old = m.cell;
      if (old >= 0) {
        auto& list = cells_[old].members;
        const int last = list.back();
        list[m.slot] = last;
        members_[last].slot = m.slot;
        list.pop_back();
        if (list.empty()) {
          cells_[old].weight_bound = 0;
        }
      }
      m.cell = cell;
      if (cell >= 0) {
        m.slot = static_cast<int>(cells_[cell].members.size());
        cells_[cell].members.push_back(person);
        if (stored_weights_) {
          cells_[cell].weight_bound =
            std::max(cells_[cell].weight_bound, m.weight);
        }
      } else {
        m.slot = -1;
      }
      if (old >= 0) {
        update(old, now);
      }
      if (cell >= 0) {
        update(cell, now);
      }
      if (on_move_) {
        on_move_(person, old, cell);
      }
    }
    boundaries_.set(person, next);
  }
  double next_boundary() { return boundaries_.top(); }
  // Move the unit with the earliest boundary to its next cell, using the
  // schedule and clock origin recorded at its last refresh.
  void advance_boundary() {
    if (++updates_ % 1024 == 0) {
      checkUserInterrupt();
    }
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
    Cell& x = cells_[c];
    if (x.members.empty()) {
      stop("A rate cell fired without members.");
    }
    const double uniform = R::runif(0, 1);
    size_t pick = uniform_index(uniform, x.members.size());
    if (stored_weights_) {
      // Stored weights: rejection sampling against the cell's weight bound,
      // exact for any bound at or above the largest weight. A member at the
      // bound is accepted without a further draw, so equal weights keep the
      // unweighted pick and its random-number use.
      // The bound only falls when the cell empties, so it can go stale after
      // a high-weight member leaves. After repeated rejections it is reset to
      // the cell's exact maximum; any bound at or above it keeps the sampling
      // exact, and with the exact maximum a pick takes at most about n tries
      // in expectation.
      for (unsigned long rejected = 0;; ++rejected) {
        const double w = members_[x.members[pick]].weight;
        if (w >= x.weight_bound || R::runif(0, 1) * x.weight_bound < w) {
          break;
        }
        if (rejected % 32 == 31) {
          double top = 0;
          for (int unit : x.members) {
            top = std::max(top, members_[unit].weight);
          }
          x.weight_bound = top;
          checkUserInterrupt();
        }
        pick = uniform_index(R::runif(0, 1), x.members.size());
      }
    } else if (chooser_) {
      pick = chooser_(c, x.members, uniform);
    }
    clock_.fired(now);
    clock_.arm(intensity_, draw_threshold_);
    return x.members[pick];
  }
  void set_chooser(
    std::function<size_t(int, const std::vector<int>&, double)> f
  ) {
    chooser_ = std::move(f);
  }
  void set_move_hook(std::function<void(int, int, int)> f) {
    on_move_ = std::move(f);
  }
  // Calendar rate changes: true if new schedules keep this process's cells
  // (same strata and nonzero band edges), so rates can change in place.
  bool same_layout(const std::vector<std::vector<Band>>& schedules) const {
    if (schedules.size() != strata_.size()) {
      return false;
    }
    for (size_t k = 0; k < schedules.size(); ++k) {
      size_t i = 0;
      for (const Band& band : schedules[k]) {
        if (!(band.rate > 0)) {
          continue;
        }
        if (i >= strata_[k].size() || band.start != strata_[k][i].start ||
            band.end != strata_[k][i].end) {
          return false;
        }
        ++i;
      }
      if (i != strata_[k].size()) {
        return false;
      }
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
        if (!(band.rate > 0)) {
          continue;
        }
        const int c = first_[k] + static_cast<int>(i);
        strata_[k][i].rate = band.rate;
        cells_[c].rate = band.rate;
        intensity_.set(
          static_cast<size_t>(c),
          band.rate * static_cast<double>(cells_[c].members.size())
        );
        ++i;
      }
    }
    clock_.arm(intensity_, draw_threshold_);
  }
  // Keep a weight per member; events choose members in proportion to it.
  void enable_stored_weights() { stored_weights_ = true; }
  void set_weight(int unit, double weight) {
    if (unit >= static_cast<int>(members_.size())) {
      members_.resize(unit + 1);
    }
    Member& m = members_[unit];
    if (m.weight == weight) {
      return;
    }
    m.weight = weight;
    if (m.cell >= 0) {
      cells_[m.cell].weight_bound =
        std::max(cells_[m.cell].weight_bound, weight);
    }
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
      if (m.cell < 0) {
        continue;
      }
      ++counts[m.cell];
      if (cells_[m.cell].members[m.slot] != static_cast<int>(p)) {
        return false;
      }
    }
    if (stored_weights_) {
      // Stored weights: each cell's bound covers its members' weights.
      for (size_t c = 0; c < cells_.size(); ++c) {
        for (int unit : cells_[c].members) {
          if (members_[unit].weight > cells_[c].weight_bound) {
            return false;
          }
        }
      }
    }
    for (size_t c = 0; c < cells_.size(); ++c) {
      if (counts[c] != static_cast<int>(cells_[c].members.size())) {
        return false;
      }
      if (intensity_.value(c) !=
          cells_[c].rate * static_cast<double>(counts[c])) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace hammel

#endif  // HAMMEL_RATE_CELLS_H
