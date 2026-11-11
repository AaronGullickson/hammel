// Helpers shared by the history-observation functions (demography.cpp and
// observation.cpp): public-ID lookup, interval bins and indexed group
// histories.
#ifndef HAMMEL_GROUP_HISTORY_H
#define HAMMEL_GROUP_HISTORY_H

#include <Rcpp.h>

#include <algorithm>
#include <map>
#include <vector>

namespace hammel {

// Positions of public person IDs in their table (the last row wins).
inline std::map<double, int> id_index(const Rcpp::NumericVector& ids) {
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) {
    index[ids[i]] = i;
  }
  return index;
}

// The half-open interval of `breaks` containing a value, or -1 outside them.
inline int break_bin(double value, const Rcpp::NumericVector& breaks) {
  if (!R_finite(value) || value < breaks[0] ||
      value >= breaks[breaks.size() - 1]) {
    return -1;
  }
  return std::upper_bound(breaks.begin(), breaks.end(), value) -
         breaks.begin() - 1;
}

// Private observation adapter: labels are indexed in R, with one row per
// change. The public population retains the original transition table.
class GroupHistory {
 public:
  struct Change {
    double time;
    int from, to;
  };

 private:
  std::vector<std::vector<Change>> changes;
  Rcpp::IntegerVector final_group;

 public:
  GroupHistory(Rcpp::DataFrame table, Rcpp::IntegerVector groups)
      : changes(table.nrows()), final_group(groups) {
    if (!table.hasAttribute(".hammel_group_history")) {
      return;
    }
    Rcpp::List history = table.attr(".hammel_group_history");
    Rcpp::IntegerVector person = history["row"], from = history["from"],
                        to = history["to"];
    Rcpp::NumericVector time = history["time"];
    for (int j = 0; j < person.size(); ++j) {
      changes[person[j]].push_back({time[j], from[j], to[j]});
    }
    for (auto& h : changes) {
      std::stable_sort(
        h.begin(), h.end(),
        [](const Change& a, const Change& b) { return a.time < b.time; }
      );
    }
  }
  int at(int row, double time, bool left = false) const {
    const auto& h = changes[row];
    if (h.empty()) {
      return final_group[row];
    }
    auto next =
      std::lower_bound(h.begin(), h.end(), time, [](const Change& c, double t) {
        return c.time < t;
      });
    // A combined history (several dimensions) can hold several changes at
    // one instant; the state after a time follows all of them.
    if (!left) {
      while (next != h.end() && next->time == time) {
        ++next;
      }
    }
    return next == h.begin() ? h.front().from : std::prev(next)->to;
  }
  const std::vector<Change>& of(int row) const { return changes[row]; }
  void append_times(int row, std::vector<double>& times) const {
    for (const auto& change : changes[row]) {
      times.push_back(change.time);
    }
  }
  template <class F>
  void split(int row, double lo, double hi, F add) const {
    if (hi <= lo) {
      return;
    }
    const auto& h = changes[row];
    auto next =
      std::upper_bound(h.begin(), h.end(), lo, [](double t, const Change& c) {
        return t < c.time;
      });
    int group = at(row, lo);
    while (next != h.end() && next->time < hi) {
      add(group, lo, next->time);
      lo = next->time;
      group = next->to;
      ++next;
    }
    add(group, lo, hi);
  }
};

}  // namespace hammel

#endif  // HAMMEL_GROUP_HISTORY_H
