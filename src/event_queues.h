// Scheduling primitives shared by the rate processes: the boundary queue,
// intensity sum trees, the process clock and weighted member selection.
#ifndef HAMMEL_EVENT_QUEUES_H
#define HAMMEL_EVENT_QUEUES_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

namespace hammel {

using namespace Rcpp;

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
// Entries are lazily deleted: one is current only while it matches
// next_[person].
class BoundaryQueue {
  static constexpr int kBuckets = 365, kSpan = 16, kRing = kBuckets * kSpan;
  using Entry = std::pair<double, int>;
  struct Later {
    bool operator()(const Entry& a, const Entry& b) const { return a > b; }
  };
  std::vector<std::vector<Entry>> ring_ =
    std::vector<std::vector<Entry>>(kRing);
  std::vector<Entry> overflow_;  // min-heap of entries beyond the ring
  std::vector<double> next_;
  size_t live_ = 0, in_ring_ = 0, pos_ = 0;
  long long cur_ = 0;  // absolute index of the bucket being drained
  bool started_ = false, sorted_ = false;
  static long long bucket(double t) {
    return static_cast<long long>(std::floor(t * kBuckets));
  }
  std::vector<Entry>& slot(long long b) {
    return ring_[static_cast<size_t>(((b % kRing) + kRing) % kRing)];
  }
  void push(const Entry& e) {
    long long b = bucket(e.first);
    if (!started_) {
      cur_ = b;
      started_ = true;
    }
    if (b < cur_) {
      b = cur_;  // at or after the current time; roundoff guard
    }
    if (b - cur_ >= kRing) {
      overflow_.push_back(e);
      std::push_heap(overflow_.begin(), overflow_.end(), Later());
      return;
    }
    auto& v = slot(b);
    if (b == cur_ && sorted_) {
      v.insert(
        std::upper_bound(v.begin() + static_cast<long>(pos_), v.end(), e), e
      );
    } else {
      v.push_back(e);
    }
    ++in_ring_;
  }
  // Position the cursor on the earliest current entry; false if none remain.
  bool settle() {
    if (live_ == 0) {
      return false;
    }
    for (;;) {
      if (!sorted_ && in_ring_ == 0 &&
          !overflow_.empty()) {  // jump over empty buckets
        cur_ = std::max(cur_, bucket(overflow_.front().first));
      }
      auto& v = slot(cur_);
      if (!sorted_) {
        while (!overflow_.empty() && bucket(overflow_.front().first) <= cur_) {
          v.push_back(overflow_.front());
          ++in_ring_;
          std::pop_heap(overflow_.begin(), overflow_.end(), Later());
          overflow_.pop_back();
        }
        std::sort(v.begin(), v.end());
        sorted_ = true;
        pos_ = 0;
      }
      while (pos_ < v.size()) {
        const Entry& e = v[pos_];
        if (e.first == next_[e.second]) {
          return true;
        }
        ++pos_;
        --in_ring_;
      }
      v.clear();
      sorted_ = false;
      ++cur_;
    }
  }

 public:
  void set(int person, double time) {
    if (person >= static_cast<int>(next_.size())) {
      next_.resize(person + 1, R_PosInf);
    }
    if (next_[person] == time) {
      return;
    }
    if (R_finite(next_[person])) {
      --live_;
    }
    next_[person] = time;
    if (R_finite(time)) {
      ++live_;
      push({time, person});
    }
  }
  void clear(int person) {
    if (person < static_cast<int>(next_.size())) {
      set(person, R_PosInf);
    }
  }
  double top() {
    if (!settle()) {
      return R_PosInf;
    }
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
};

// A binary sum tree over growable leaves: leaf i holds a nonnegative value
// and internal nodes hold the sums of their children. Sums are recomputed
// from children along each updated path, so they do not drift, and removing
// a dominant value restores small ones accurately. Capacity doubles as
// leaves are added.
class SumTree {
 public:
  // Set leaf i, growing the tree to hold it.
  void set(size_t i, double value) {
    if (i >= capacity_) {
      grow(i);
    }
    size_t node = capacity_ + i;
    tree_[node] = value;
    while (node > 1) {
      node /= 2;
      tree_[node] = tree_[2 * node] + tree_[2 * node + 1];
    }
  }
  double total() const { return tree_.empty() ? 0 : tree_[1]; }
  double value(size_t i) const {
    return i < capacity_ ? tree_[capacity_ + i] : 0;
  }

 protected:
  size_t capacity_ = 0;
  std::vector<double> tree_;

 private:
  void grow(size_t i) {
    size_t next = capacity_ ? capacity_ : 1;
    while (next <= i) {
      next *= 2;
    }
    std::vector<double> grown(2 * next, 0);
    for (size_t k = 0; k < capacity_; ++k) {
      grown[next + k] = tree_[capacity_ + k];
    }
    for (size_t k = next - 1; k > 0; --k) {
      grown[k] = grown[2 * k] + grown[2 * k + 1];
    }
    capacity_ = next;
    tree_.swap(grown);
  }
};

// Event intensities of a process's cells (leaf c is cell c).
class IntensityTree : public SumTree {
 public:
  // The cell holding position uniform * total; only cells with positive
  // intensity.
  size_t choose(double uniform) const {
    double value = std::min(uniform * total(), std::nextafter(total(), 0.0));
    size_t node = 1;
    while (node < capacity_) {
      const double left = tree_[2 * node], right = tree_[2 * node + 1];
      // Roundoff guard: never descend into a branch with no intensity.
      if ((value < left && left > 0) || !(right > 0)) {
        node *= 2;
      } else {
        value = std::max(0.0, value - left);
        node = 2 * node + 1;
      }
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
    if (has_threshold_ && now > last_) {
      threshold_ = std::max(0.0, threshold_ - tree.total() * (now - last_));
    }
    last_ = now;
  }
  // After intensity changes: ensure a threshold exists while events can occur.
  void arm(const IntensityTree& tree, const std::function<double()>& draw) {
    if (!has_threshold_ && tree.total() > 0) {
      threshold_ = draw();
      has_threshold_ = true;
    }
  }
  double due(const IntensityTree& tree) const {
    return has_threshold_ && tree.total() > 0
             ? last_ + threshold_ / tree.total()
             : R_PosInf;
  }
  // The clock fired at now: its threshold is used up.
  void fired(double now) {
    has_threshold_ = false;
    last_ = now;
  }
};

// The member floor(uniform * n) of n, guarded against uniform rounding up to
// one: the unweighted pick shared by every process.
inline size_t uniform_index(double uniform, size_t n) {
  return std::min(
    n - 1, static_cast<size_t>(std::floor(uniform * static_cast<double>(n)))
  );
}

// Choose a member from log weights with one uniform draw. Equal log weights
// use the unweighted rule floor(uniform * n), so constant weights reproduce
// unweighted runs exactly (D129).
inline size_t pick_log_weights(const std::vector<double>& eta, double uniform) {
  const size_t n = eta.size();
  bool equal = true;
  double top = eta[0];
  for (double e : eta) {
    if (e != eta[0]) {
      equal = false;
    }
    top = std::max(top, e);
  }
  if (equal) {
    return uniform_index(uniform, n);
  }
  std::vector<double> cumulative(n);
  double total = 0;
  for (size_t k = 0; k < n; ++k) {
    total += std::exp(eta[k] - top);
    cumulative[k] = total;
  }
  const double target = uniform * total;
  for (size_t k = 0; k < n; ++k) {
    if (target < cumulative[k]) {
      return k;
    }
  }
  return n - 1;
}

// The selection weights of a cell's eligible women, in the same swap-remove
// slots as the cell's eligible list.
class FertilityWeightTree : public SumTree {
 public:
  void add(double value) { put(size_++, value); }
  void remove(size_t slot) {
    const size_t last = size_ - 1;
    const double moved = value(last);
    put(slot, 0);
    if (slot != last) {
      put(last, 0);
      put(slot, moved);
    }
    --size_;
  }
  // The slot holding position uniform * total.
  size_t choose(double uniform) const {
    double value = std::min(uniform * total(), std::nextafter(total(), 0.0));
    size_t node = 1;
    while (node < capacity_) {
      const double left = tree_[2 * node];
      if (value < left) {
        node *= 2;
      } else {
        value -= left;
        node = 2 * node + 1;
      }
    }
    return std::min(node - capacity_, size_ - 1);
  }

 private:
  size_t size_ = 0;
  void put(size_t slot, double value) {
    set(slot, value);
    if (!R_finite(total())) {
      stop("Eligible hetfert weight sum is not representable.");
    }
  }
};

}  // namespace hammel

#endif  // HAMMEL_EVENT_QUEUES_H
