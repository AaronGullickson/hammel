#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <vector>
#include "group_history.h"

using namespace Rcpp;

// Shared observation engine (D133). Each person has a timeline of states that
// change only at recorded events. Same-time changes are ordered by rank, which
// follows the simulation's event order; an event at rank r sees the changes
// at its own time whose rank is lower.
namespace {

int observation_bin(double value, const NumericVector& breaks) {
  if (!R_finite(value) || value < breaks[0] || value >= breaks[breaks.size() - 1])
    return -1;
  return std::upper_bound(breaks.begin(), breaks.end(), value) - breaks.begin() - 1;
}

// Change ranks. A union start that touches the end of the person's previous
// union at the same instant takes that end's rank plus one half, so births at
// that instant still see the ending union (birth precedes dissolution and
// conversion).
const double kDeathEnd = 0, kStart = 1, kBirth = 2, kDissolutionEnd = 3,
  kConversionEnd = 4, kTransition = 6;
// Event ranks.
const double kDeathEvent = 0, kBirthEvent = 2, kFormationEvent = 5, kTransitionEvent = 6;

struct State {
  int group = 0, utype = 0, previous = 0, legal = 0, parity = 0;
  double group_entry = 0, union_end = NA_REAL;
};

struct Change {
  double time, rank;
  int kind;      // 0 group, 1 union start, 2 union end, 3 birth
  int value;     // group index, union type, previous-union category
  int legal;     // legal previous category after an ending
};

struct Snapshot { double time, rank; State state; };

struct Episode {
  double start, end, id;
  int after;           // 1 dissolution, 2 other ending, 0 ongoing
  bool cohab, conversion, death_end;
};

struct Formation { double time; int kind; };  // 0 cohab, 1 marriage, 2 converted

// Inputs shared by the counters, read once from the R tables.
struct Histories {
  int n;
  NumericVector ids, born, died, mothers;
  std::vector<int> female;
  std::vector<std::vector<Snapshot>> timelines;
  std::vector<std::vector<Formation>> formations;
  std::vector<std::vector<int>> children;  // child rows, by mother
  std::vector<int> mother_row;
  std::vector<std::vector<std::pair<double, int>>> moves;  // time, destination
};

Histories build_histories(DataFrame persons, DataFrame unions, IntegerVector groups,
                          bool typed, bool exclusive, bool counts_cohab) {
  Histories h;
  h.ids = persons["person_id"]; h.born = persons["birth_time"];
  h.died = persons["death_time"]; h.mothers = persons["mother_id"];
  h.n = h.ids.size();
  CharacterVector sex = persons["sex"];
  h.female.resize(h.n);
  std::map<double, int> index;
  for (int i = 0; i < h.n; ++i) {
    index[h.ids[i]] = i;
    h.female[i] = sex[i] == "f";
  }
  h.mother_row.assign(h.n, -1);
  h.children.resize(h.n);
  for (int i = 0; i < h.n; ++i) if (!NumericVector::is_na(h.mothers[i])) {
    auto m = index.find(h.mothers[i]);
    if (m == index.end()) stop("A mother is missing from this history group.");
    if (!h.female[m->second]) stop("Birth-rate extraction requires female maternal records.");
    h.mother_row[i] = m->second;
    h.children[m->second].push_back(i);
  }
  for (auto& c : h.children)
    std::stable_sort(c.begin(), c.end(), [&](int a, int b) { return h.born[a] < h.born[b]; });

  // Union episodes for each partner.
  std::vector<std::vector<Episode>> episodes(h.n);
  if (unions.nrows()) {
    NumericVector wives = unions["wife_id"], husbands = unions["husband_id"];
    NumericVector starts = unions["start_time"], ends = unions["end_time"], uids = unions["union_id"];
    CharacterVector reasons = unions["end_reason"];
    CharacterVector types = typed ? as<CharacterVector>(unions["union_type"]) : CharacterVector();
    for (int u = 0; u < starts.size(); ++u) {
      const bool ongoing = NumericVector::is_na(ends[u]);
      const std::string reason = ongoing || CharacterVector::is_na(reasons[u]) ? "" :
        as<std::string>(reasons[u]);
      const double end = ongoing ? R_PosInf : ends[u];
      for (double id : {wives[u], husbands[u]}) {
        const int i = index.at(id);
        if (exclusive && !NumericVector::is_na(h.died[i]) && (starts[u] > h.died[i] || end > h.died[i]))
          stop("Previous-union histories cannot contain union events after a partner's death.");
        episodes[i].push_back({starts[u], end, uids[u], ongoing ? 0 : (reason == "dissolution" ? 1 : 2),
          typed && types[u] == "cohab", reason == "conversion",
          reason == "death_husband" || reason == "death_wife"});
      }
    }
  }
  GroupHistory group_history(persons, groups);
  h.timelines.resize(h.n);
  h.formations.resize(h.n);
  h.moves.resize(h.n);
  for (int i = 0; i < h.n; ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    auto& e = episodes[i];
    std::sort(e.begin(), e.end(), [](const Episode& a, const Episode& b) {
      return a.start < b.start || (a.start == b.start && a.id < b.id);
    });
    if (exclusive) {
      for (size_t u = 1; u < e.size(); ++u)
        if (e[u].start < e[u - 1].end)
          stop("Observation by union type or previous union requires nonoverlapping union histories for each person.");
    } else {
      // Untyped coverage only: merge overlapping episodes. A merged ending is
      // birth-first if any merged component ended by dissolution there.
      std::vector<Episode> merged;
      for (const auto& x : e) {
        if (!merged.empty() && x.start <= merged.back().end) {
          auto& m = merged.back();
          if (x.end > m.end) { m.end = x.end; m.after = x.after; m.death_end = x.death_end; }
          else if (x.end == m.end && x.after == 1) m.after = 1;
        } else merged.push_back(x);
      }
      e = std::move(merged);
    }
    std::vector<Change> changes;
    double previous_end = R_NegInf, previous_end_rank = 0;
    bool converted = false;
    for (const auto& x : e) {
      const double start_rank = x.start == previous_end ? previous_end_rank + 0.5 : kStart;
      changes.push_back({x.start, start_rank, 1, x.cohab ? 1 : 2, 0});
      h.formations[i].push_back({x.start, converted && x.start == previous_end ? 2 : (x.cohab ? 0 : 1)});
      if (R_finite(x.end)) {
        const double rank = x.conversion ? kConversionEnd : (x.after == 1 ? kDissolutionEnd :
          (x.death_end ? kDeathEnd : kDissolutionEnd));
        changes.push_back({x.end, rank, 2, x.after, x.conversion ? -1 : (x.cohab ? -2 : x.after)});
        previous_end = x.end;
        previous_end_rank = rank;
        converted = x.conversion;
      }
    }
    for (int child : h.children[i]) changes.push_back({h.born[child], kBirth, 3, 0, 0});
    std::vector<double> times;
    group_history.append_times(i, times);
    for (double t : times) {
      const int to = group_history.at(i, t);
      changes.push_back({t, kTransition, 0, to, 0});
      h.moves[i].push_back({t, to});
    }
    std::stable_sort(changes.begin(), changes.end(), [](const Change& a, const Change& b) {
      return a.time < b.time || (a.time == b.time && a.rank < b.rank);
    });
    State s;
    s.group = times.empty() ? groups[i] : group_history.at(i, h.born[i], true);
    s.group_entry = h.born[i];
    auto& line = h.timelines[i];
    line.push_back({R_NegInf, -1, s});
    for (const auto& c : changes) {
      switch (c.kind) {
        case 0: s.group = c.value; s.group_entry = c.time; break;
        case 1: s.utype = typed ? c.value : 1; break;
        case 2:
          s.utype = 0;
          s.union_end = c.time;
          // Conversion is not an ending for previous-union purposes.
          if (c.legal != -1) {
            const bool cohab = c.legal == -2;
            if (!cohab) s.legal = c.value;
            s.previous = cohab && !counts_cohab ? s.legal : c.value;
          }
          break;
        case 3: ++s.parity; break;
      }
      line.push_back({c.time, c.rank, s});
    }
  }
  return h;
}

// State seen by an event at (time, rank).
const State& state_at(const std::vector<Snapshot>& line, double time, double rank) {
  auto next = std::lower_bound(line.begin() + 1, line.end(), std::make_pair(time, rank),
    [](const Snapshot& s, const std::pair<double, double>& key) {
      return s.time < key.first || (s.time == key.first && s.rank < key.second);
    });
  return std::prev(next)->state;
}

struct Layout {
  int nt, na, ns, np, ng, nk, nr;
  bool sex, parity, group, utype, previous;
  IntegerVector previous_map;
  R_xlen_t cells() const { return static_cast<R_xlen_t>(nt) * na * ns * np * ng * nk * nr; }
  R_xlen_t cell(int t, int a, int s, int p, const State& x) const {
    const int g = group ? x.group : 0, k = utype ? x.utype : 0,
      r = previous ? previous_map[x.previous] : 0;
    return t + static_cast<R_xlen_t>(nt) * (a + static_cast<R_xlen_t>(na) * (s +
      static_cast<R_xlen_t>(ns) * (p + static_cast<R_xlen_t>(np) * (g +
      static_cast<R_xlen_t>(ng) * (k + static_cast<R_xlen_t>(nk) * r)))));
  }
};

}  // namespace

// Person-level Lexis counts. Exposure and events are split at every state
// change. `axis`: 0 age, 1 years since group entry (birth before any move),
// 2 years since the most recent union ending (persons without one are
// skipped). `utype_filter`: 0 all time, 1 unpartnered only, 2 unpartnered or
// cohabiting. `previous_mask` limits exposure and formations to the selected
// previous-union categories.
// [[Rcpp::export]]
List person_lexis_cpp(DataFrame persons, DataFrame unions, IntegerVector groups,
    int group_count, NumericVector time_breaks, NumericVector axis_breaks,
    NumericVector parity_breaks, List options) {
  const bool typed = options["typed"], exclusive = options["exclusive"],
    counts_cohab = options["counts_cohab"];
  const int axis = options["axis"], utype_filter = options["utype_filter"];
  LogicalVector previous_mask = options["previous_mask"];
  Layout L;
  L.sex = options["sex"]; L.group = options["group"]; L.utype = options["union_type"];
  L.previous = options["previous_union"]; L.parity = parity_breaks.size() > 1;
  L.previous_map = options["previous_map"];
  L.nt = time_breaks.size() - 1; L.na = axis_breaks.size() - 1;
  L.ns = L.sex ? 2 : 1; L.np = L.parity ? parity_breaks.size() - 1 : 1;
  L.ng = L.group ? group_count : 1; L.nk = L.utype ? (typed ? 3 : 2) : 1;
  L.nr = L.previous ? max(L.previous_map) + 1 : 1;
  const double horizon = time_breaks[L.nt];
  Histories h = build_histories(persons, unions, groups, typed, exclusive, counts_cohab);
  const R_xlen_t cells = L.cells();
  NumericVector exposure(cells), deaths(cells), births(cells), formations(cells * 3),
    transitions(cells * group_count), unknown(L.nt);
  auto parity_bin = [&](int parity) {
    return L.parity ? observation_bin(static_cast<double>(parity), parity_breaks) : 0;
  };
  auto origin_of = [&](int i, const State& x) {
    return axis == 0 ? h.born[i] : (axis == 1 ? x.group_entry : x.union_end);
  };
  auto included = [&](const State& x) {
    if (utype_filter == 1 && x.utype != 0) return false;
    if (utype_filter == 2 && x.utype == (typed ? 2 : 1)) return false;
    if (!previous_mask[x.previous]) return false;
    if (axis == 2 && !R_finite(x.union_end)) return false;
    return true;
  };
  for (int i = 0; i < h.n; ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    const int s = L.sex && h.female[i] ? 1 : 0;
    const double life_end = NumericVector::is_na(h.died[i]) ? horizon : h.died[i];
    const auto& line = h.timelines[i];
    // Exposure: one piece per stretch of constant state.
    for (size_t j = 0; j < line.size(); ++j) {
      const double lo = std::max(j ? line[j].time : h.born[i], h.born[i]);
      const double hi = std::min(j + 1 < line.size() ? line[j + 1].time : life_end, life_end);
      if (hi <= lo) continue;
      const State& x = line[j].state;
      if (!included(x)) continue;
      const int p = parity_bin(x.parity);
      if (p < 0) continue;
      const double origin = origin_of(i, x);
      const double from = std::max(lo, time_breaks[0]), to = std::min(hi, horizon);
      if (to <= from) continue;
      int t = observation_bin(from, time_breaks);
      for (; t >= 0 && t < L.nt && time_breaks[t] < to; ++t) {
        const double a_lo = std::max(from, time_breaks[t]), a_hi = std::min(to, time_breaks[t + 1]);
        if (a_hi <= a_lo) continue;
        int a = observation_bin(a_lo - origin, axis_breaks);
        if (a < 0) {
          // Before the first axis band: start at the first band.
          if (a_lo - origin < axis_breaks[0]) a = 0; else continue;
        }
        for (; a < L.na; ++a) {
          const double enter = std::max(a_lo, origin + axis_breaks[a]);
          const double exit = std::min(a_hi, origin + axis_breaks[a + 1]);
          if (exit > enter) exposure[L.cell(t, a, s, p, x)] += exit - enter;
          if (origin + axis_breaks[a + 1] >= a_hi) break;
        }
      }
    }
    // Death.
    if (!NumericVector::is_na(h.died[i])) {
      const State& x = state_at(line, h.died[i], kDeathEvent);
      const int t = observation_bin(h.died[i], time_breaks);
      const int a = observation_bin(h.died[i] - origin_of(i, x), axis_breaks);
      const int p = parity_bin(x.parity);
      if (t >= 0 && a >= 0 && p >= 0 && included(x)) deaths[L.cell(t, a, s, p, x)] += 1;
    }
    // Births to this mother, each at the parity before it.
    const auto& kids = h.children[i];
    for (size_t j = 0; j < kids.size(); ++j) {
      const double when = h.born[kids[j]];
      State x = state_at(line, when, kBirthEvent);
      x.parity = static_cast<int>(j);
      const int t = observation_bin(when, time_breaks);
      const int a = observation_bin(when - origin_of(i, x), axis_breaks);
      const int p = parity_bin(x.parity);
      if (t >= 0 && a >= 0 && p >= 0 && included(x)) births[L.cell(t, a, L.sex ? 1 : 0, p, x)] += 1;
    }
    // Union formations, under the state before them.
    for (const auto& f : h.formations[i]) {
      const State& x = state_at(line, f.time, kFormationEvent);
      if (!previous_mask[x.previous]) continue;
      if (axis == 2 && !R_finite(x.union_end)) continue;
      const int t = observation_bin(f.time, time_breaks);
      const int a = observation_bin(f.time - origin_of(i, x), axis_breaks);
      const int p = parity_bin(x.parity);
      if (t >= 0 && a >= 0 && p >= 0) {
        State y = x;
        y.utype = 0;
        formations[L.cell(t, a, s, p, y) + cells * f.kind] += 1;
      }
    }
    // Group transitions, from the origin group.
    for (const auto& m : h.moves[i]) {
      const State& x = state_at(line, m.first, kTransitionEvent);
      const int t = observation_bin(m.first, time_breaks);
      const int a = observation_bin(m.first - origin_of(i, x), axis_breaks);
      const int p = parity_bin(x.parity);
      if (t >= 0 && a >= 0 && p >= 0 && included(x))
        transitions[L.cell(t, a, s, p, x) + cells * m.second] += 1;
    }
  }
  // Births with unknown mothers, by period.
  for (int i = 0; i < h.n; ++i) if (h.mother_row[i] < 0) {
    const int t = observation_bin(h.born[i], time_breaks);
    if (t >= 0) unknown[t] += 1;
  }
  return List::create(_["person_years"] = exposure, _["deaths"] = deaths,
    _["births"] = births, _["formations"] = formations, _["transitions"] = transitions,
    _["unassigned_births"] = unknown);
}

// Living persons at time points, by state after all events at each time.
// [[Rcpp::export]]
NumericVector person_snapshot_cpp(DataFrame persons, DataFrame unions, IntegerVector groups,
    int group_count, NumericVector times, List options) {
  const bool typed = options["typed"], exclusive = options["exclusive"],
    counts_cohab = options["counts_cohab"];
  Layout L;
  L.sex = options["sex"]; L.group = options["group"]; L.utype = options["union_type"];
  L.previous = options["previous_union"]; L.parity = false;
  L.previous_map = options["previous_map"];
  L.nt = times.size(); L.na = 1; L.np = 1;
  L.ns = L.sex ? 2 : 1; L.ng = L.group ? group_count : 1; L.nk = L.utype ? (typed ? 3 : 2) : 1;
  L.nr = L.previous ? max(L.previous_map) + 1 : 1;
  Histories h = build_histories(persons, unions, groups, typed, exclusive, counts_cohab);
  NumericVector counts(L.cells());
  for (int i = 0; i < h.n; ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    const int s = L.sex && h.female[i] ? 1 : 0;
    for (int t = 0; t < L.nt; ++t) {
      const double when = times[t];
      if (h.born[i] > when || (!NumericVector::is_na(h.died[i]) && h.died[i] <= when)) continue;
      counts[L.cell(t, 0, s, 0, state_at(h.timelines[i], when, R_PosInf))] += 1;
    }
  }
  return counts;
}
