// Native counters for union-level demographic measures and for union and
// birth histories used by the R analysis functions.
#include <Rcpp.h>

#include <algorithm>
#include <map>

#include "group_history.h"

using namespace Rcpp;
using hammel::break_bin;
using hammel::GroupHistory;
using hammel::id_index;

// Union episodes, rather than persons: each record contributes exposure once.
// Death ends its exposure but only non-death dissolution enters the numerator.
// [[Rcpp::export]]
List dissolution_counts_cpp(
  DataFrame unions,
  NumericVector time_breaks,
  NumericVector dur_breaks,
  IntegerVector group_indices,
  int group_count,
  Nullable<NumericVector> origins = R_NilValue,
  std::string event_reason = "dissolution"
) {
  GroupHistory group_history(unions, group_indices);
  NumericVector starts = unions["start_time"], ends = unions["end_time"];
  const NumericVector axis_origins =
    origins.isNotNull() ? NumericVector(origins) : starts;
  CharacterVector reasons = unions["end_reason"];
  const int nt = time_breaks.size() - 1, nd = dur_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * nd;
  NumericVector exposure(stride * group_count), events(stride * group_count);
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) {
      checkUserInterrupt();
    }
    const double end =
      NumericVector::is_na(ends[u]) ? time_breaks[nt] : ends[u];
    for (int t = 0; t < nt; ++t) {
      const double lo = std::max(starts[u], time_breaks[t]);
      const double hi = std::min(end, time_breaks[t + 1]);
      if (hi <= lo) {
        continue;
      }
      for (int d = 0; d < nd; ++d) {
        const double enter = std::max(lo, axis_origins[u] + dur_breaks[d]);
        const double exit = std::min(hi, axis_origins[u] + dur_breaks[d + 1]);
        group_history.split(u, enter, exit, [&](int g, double lo, double hi) {
          exposure[stride * g + t + nt * d] += hi - lo;
        });
      }
    }
    if (NumericVector::is_na(ends[u]) || CharacterVector::is_na(reasons[u]) ||
        as<std::string>(reasons[u]) != event_reason) {
      continue;
    }
    const int t = break_bin(ends[u], time_breaks);
    const int d = break_bin(ends[u] - axis_origins[u], dur_breaks);
    if (t >= 0 && d >= 0) {
      events[stride * group_history.at(u, ends[u], true) + t + nt * d] += 1;
    }
  }
  return List::create(_["union_years"] = exposure, _["dissolutions"] = events);
}

// Exclusive histories for previous-union state. Sort by time and union ID so
// input row order never changes status, including touching/zero-length
// episodes.
struct MarriageEpisode {
  double start, end, id;
  int after;
  bool cohab, conversion;
};
std::vector<std::vector<MarriageEpisode>> marriage_histories(
  DataFrame persons,
  DataFrame unions,
  const std::string& context = "Previous-union stratification"
) {
  NumericVector ids = persons["person_id"], died = persons["death_time"];
  NumericVector wives = unions["wife_id"], husbands = unions["husband_id"];
  NumericVector starts = unions["start_time"], ends = unions["end_time"],
                uids = unions["union_id"];
  CharacterVector reasons = unions["end_reason"];
  const bool typed = unions.containsElementNamed("union_type");
  CharacterVector types;
  if (typed) {
    types = unions["union_type"];
  }
  const std::map<double, int> index = id_index(ids);
  std::vector<std::vector<MarriageEpisode>> histories(ids.size());
  for (int u = 0; u < starts.size(); ++u) {
    const bool ongoing = NumericVector::is_na(ends[u]);
    const double end = ongoing ? R_PosInf : ends[u];
    const int after =
      !ongoing && as<std::string>(reasons[u]) == "dissolution" ? 1 : 2;
    for (double id : {wives[u], husbands[u]}) {
      const int i = index.at(id);
      if (!NumericVector::is_na(died[i]) &&
          (starts[u] > died[i] || end > died[i])) {
        stop(
          "Previous-union histories cannot contain union events after a "
          "partner's death."
        );
      }
      histories[i].push_back(
        {starts[u], end, uids[u], after,
         typed && as<std::string>(types[u]) == "cohab",
         !ongoing && as<std::string>(reasons[u]) == "conversion"}
      );
    }
  }
  for (auto& history : histories) {
    std::sort(
      history.begin(), history.end(),
      [](const MarriageEpisode& a, const MarriageEpisode& b) {
        return a.start < b.start || (a.start == b.start && a.id < b.id);
      }
    );
    for (size_t u = 1; u < history.size(); ++u) {
      if (history[u].start < history[u - 1].end) {
        stop(
          context + " requires nonoverlapping union histories for each person."
        );
      }
    }
  }
  return histories;
}

// Execution-only cache: complete supplied history determines the latest ending.
// [[Rcpp::export]]
List previous_union_info_cpp(
  DataFrame persons,
  DataFrame unions,
  bool counts_cohab = false
) {
  const auto histories = marriage_histories(persons, unions);
  IntegerVector state(histories.size(), 0), legal(histories.size(), 0);
  NumericVector ended(histories.size(), NA_REAL);
  for (size_t i = 0; i < histories.size(); ++i) {
    for (const auto& episode : histories[i]) {
      if (R_finite(episode.end) && !episode.conversion) {
        if (!episode.cohab) {
          legal[i] = episode.after;
        }
        state[i] = episode.cohab && !counts_cohab ? legal[i] : episode.after;
        ended[i] = episode.end;
      }
    }
  }
  return List::create(
    _["previous_union"] = state, _["end_time"] = ended,
    _["legal_previous"] = legal
  );
}

// Existing internal snapshot interface for category-only validation.
// [[Rcpp::export]]
IntegerVector previous_union_state_cpp(DataFrame persons, DataFrame unions) {
  return as<IntegerVector>(
    previous_union_info_cpp(persons, unions)["previous_union"]
  );
}

// Row selection for parent-enriched father-choice observations. Merge union
// coverage so duplicate/overlapping episodes never duplicate births.
// [[Rcpp::export]]
IntegerVector unpartnered_birth_rows_cpp(DataFrame persons, DataFrame unions) {
  NumericVector ids = persons["person_id"], born = persons["birth_time"],
                mothers = persons["mother_id"], fathers = persons["father_id"];
  NumericVector wives = unions["wife_id"], starts = unions["start_time"],
                ends = unions["end_time"];
  CharacterVector reasons = unions["end_reason"];
  const std::map<double, int> index = id_index(ids);
  std::vector<std::vector<std::pair<double, double>>> coverage(ids.size());
  std::vector<std::vector<double>> birth_first_ends(ids.size());
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) {
      checkUserInterrupt();
    }
    const int mother = index.at(wives[u]);
    coverage[mother].push_back(
      {starts[u], NumericVector::is_na(ends[u]) ? R_PosInf : ends[u]}
    );
    if (!NumericVector::is_na(ends[u]) && !CharacterVector::is_na(reasons[u])) {
      const std::string reason = as<std::string>(reasons[u]);
      if (reason == "dissolution" || reason == "conversion") {
        birth_first_ends[mother].push_back(ends[u]);
      }
    }
  }
  for (size_t i = 0; i < coverage.size(); ++i) {
    auto& intervals = coverage[i];
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<double, double>> merged;
    for (const auto& interval : intervals) {
      if (merged.empty() || interval.first > merged.back().second) {
        merged.push_back(interval);
      } else {
        merged.back().second = std::max(merged.back().second, interval.second);
      }
    }
    intervals = std::move(merged);
    std::sort(birth_first_ends[i].begin(), birth_first_ends[i].end());
  }
  std::vector<int> rows;
  for (int child = 0; child < ids.size(); ++child) {
    if (child % 1024 == 0) {
      checkUserInterrupt();
    }
    if (NumericVector::is_na(mothers[child]) ||
        NumericVector::is_na(fathers[child])) {
      continue;
    }
    const int mother = index.at(mothers[child]);
    const auto& intervals = coverage[mother];
    const auto after = std::upper_bound(
      intervals.begin(), intervals.end(), born[child],
      [](double time, const std::pair<double, double>& interval) {
        return time < interval.first;
      }
    );
    const auto& ties = birth_first_ends[mother];
    const bool partnered =
      (after != intervals.begin() && born[child] < std::prev(after)->second) ||
      std::binary_search(ties.begin(), ties.end(), born[child]);
    if (!partnered) {
      rows.push_back(child + 1);
    }
  }
  return wrap(rows);
}

// Partner identity at formation, before a same-time group transition. Reuse
// indexed histories to avoid repeated R scans and named-list lookups per union.
// [[Rcpp::export]]
IntegerMatrix union_start_groups_cpp(
  DataFrame persons,
  DataFrame unions,
  IntegerVector group_indices
) {
  GroupHistory history(persons, group_indices);
  NumericVector ids = persons["person_id"], wives = unions["wife_id"],
                husbands = unions["husband_id"], starts = unions["start_time"];
  const std::map<double, int> index = id_index(ids);
  IntegerMatrix groups(unions.nrows(), 2);
  for (int u = 0; u < unions.nrows(); ++u) {
    if (u % 1024 == 0) {
      checkUserInterrupt();
    }
    groups(u, 0) = history.at(index.at(wives[u]), starts[u], true);
    groups(u, 1) = history.at(index.at(husbands[u]), starts[u], true);
  }
  return groups;
}

// Construct union-level indexed group histories in one native pass. Each union
// visits only its partners' changes, never a named R list of all movers.
// `mode`: 0 wife's group, 1 husband's group, 2 both (wife + G * husband).
// [[Rcpp::export]]
List union_group_history_cpp(
  DataFrame persons,
  DataFrame unions,
  IntegerVector group_indices,
  int mode,
  bool typed,
  int identity_count,
  int group_count
) {
  GroupHistory history(persons, group_indices);
  NumericVector ids = persons["person_id"], wives = unions["wife_id"],
                husbands = unions["husband_id"];
  CharacterVector types =
    typed ? as<CharacterVector>(unions["union_type"]) : CharacterVector();
  const std::map<double, int> index = id_index(ids);
  std::vector<int> rows, from, to;
  std::vector<double> times;
  for (int u = 0; u < unions.nrows(); ++u) {
    if (u % 1024 == 0) {
      checkUserInterrupt();
    }
    const int wife = index.at(wives[u]), husband = index.at(husbands[u]);
    std::vector<double> changes;
    if (mode != 1) {
      history.append_times(wife, changes);
    }
    if (mode != 0) {
      history.append_times(husband, changes);
    }
    std::sort(changes.begin(), changes.end());
    changes.erase(std::unique(changes.begin(), changes.end()), changes.end());
    const int offset = typed && types[u] == "marriage" ? group_count : 0;
    auto combine = [&](double time, bool left) {
      const int w = history.at(wife, time, left),
                h = history.at(husband, time, left);
      return mode == 0 ? w : (mode == 1 ? h : w + identity_count * h);
    };
    for (double time : changes) {
      int before = combine(time, true), after = combine(time, false);
      rows.push_back(u);
      times.push_back(time);
      from.push_back(before + offset);
      to.push_back(after + offset);
    }
  }
  return List::create(
    _["row"] = wrap(rows), _["time"] = wrap(times), _["from"] = wrap(from),
    _["to"] = wrap(to)
  );
}
