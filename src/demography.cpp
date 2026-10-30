#include <Rcpp.h>
#include <algorithm>
#include <map>
#include "group_history.h"

using namespace Rcpp;

// Find the half-open interval containing a value, or -1 outside the breaks.
int demographic_bin(double value, const NumericVector& breaks) {
  if (!R_finite(value) || value < breaks[0] || value >= breaks[breaks.size() - 1])
    return -1;
  return std::upper_bound(breaks.begin(), breaks.end(), value) - breaks.begin() - 1;
}

// [[Rcpp::export]]
List demographic_counts_cpp(DataFrame persons, NumericVector time_breaks,
                             NumericVector age_breaks, IntegerVector group_indices,
                             int group_count) {
  GroupHistory group_history(persons, group_indices);
  NumericVector ids = persons["person_id"], born = persons["birth_time"];
  NumericVector died = persons["death_time"], mothers = persons["mother_id"];
  CharacterVector sex = persons["sex"];
  int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * na * 2;
  const R_xlen_t cells = stride * group_count;
  NumericVector exposure(cells), deaths(cells), births(cells), unknown(nt);
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  for (int i = 0; i < ids.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
    double end = NumericVector::is_na(died[i]) ? time_breaks[nt] : died[i];
    // Clip each life to the period, then split it at the actual age boundaries.
    for (int t = 0; t < nt; ++t) {
      double lo = std::max(born[i], time_breaks[t]);
      double hi = std::min(end, time_breaks[t + 1]);
      if (hi <= lo) continue;
      for (int a = 0; a < na; ++a) {
        double enter = std::max(lo, born[i] + age_breaks[a]);
        double exit = std::min(hi, born[i] + age_breaks[a + 1]);
        group_history.split(i, enter, exit, [&](int g, double lo, double hi) {
          exposure[t + nt * (a + na * s) + stride * g] += hi - lo;
        });
      }
    }
    int t = demographic_bin(died[i], time_breaks);
    int a = demographic_bin(died[i] - born[i], age_breaks);
    if (t >= 0 && a >= 0) deaths[t + nt * (a + na * s) + stride * group_history.at(i, died[i], true)] += 1;

    t = demographic_bin(born[i], time_breaks);
    if (t < 0) continue;
    if (NumericVector::is_na(mothers[i])) {
      unknown[t] += 1;
      continue;
    }
    auto mother = index.find(mothers[i]);
    if (mother == index.end()) stop("A mother is missing from this history group.");
    if (as<std::string>(sex[mother->second]) != "f")
      stop("Birth-rate extraction requires female maternal records.");
    a = demographic_bin(born[i] - born[mother->second], age_breaks);
    if (a >= 0) births[t + nt * (a + na) + stride * group_history.at(mother->second, born[i], true)] += 1;
  }
  return List::create(_["person_years"] = exposure, _["deaths"] = deaths,
                       _["births"] = births, _["unassigned_births"] = unknown);
}

// Optional union summaries use union-history coverage, not sums of episodes:
// overlapping records must never count the same person-time twice.
// [[Rcpp::export]]
List union_counts_cpp(DataFrame persons, DataFrame unions,
                      NumericVector time_breaks, NumericVector age_breaks, IntegerVector group_indices,
                             int group_count) {
  GroupHistory group_history(persons, group_indices);
  NumericVector ids = persons["person_id"], born = persons["birth_time"];
  NumericVector died = persons["death_time"];
  CharacterVector sex = persons["sex"];
  NumericVector wives = unions["wife_id"], husbands = unions["husband_id"];
  NumericVector starts = unions["start_time"], ends = unions["end_time"];
  CharacterVector end_reasons = unions["end_reason"];
  const int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * na * 2;
  const R_xlen_t cells = stride * group_count;
  NumericVector partnered(cells);
  std::map<double, int> index;
  std::vector<std::vector<std::pair<double, double>>> episodes(ids.size());
  // Dissolution follows a tied birth (D078). Keep these exact endpoints
  // separately: merging intervals is correct for exposure, but loses reasons.
  std::vector<std::vector<double>> dissolution_ends(ids.size());
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    if (!NumericVector::is_na(ends[u]) &&
        !CharacterVector::is_na(end_reasons[u]) &&
        as<std::string>(end_reasons[u]) == "dissolution")
      dissolution_ends[index.at(wives[u])].push_back(ends[u]);
    for (double id : {wives[u], husbands[u]}) {
      const int i = index.at(id);
      const double end = NumericVector::is_na(ends[u]) ? time_breaks[nt] : ends[u];
      episodes[i].push_back({starts[u], end});

    }
  }
  for (int i = 0; i < ids.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    auto& dissolution_times = dissolution_ends[i];
    std::sort(dissolution_times.begin(), dissolution_times.end());
    auto& intervals = episodes[i];
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<double, double>> merged;
    for (const auto& interval : intervals) {
      if (merged.empty() || interval.first > merged.back().second)
        merged.push_back(interval);
      else merged.back().second = std::max(merged.back().second, interval.second);
    }
    intervals = std::move(merged);
    const int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
    const double life_end = NumericVector::is_na(died[i]) ? time_breaks[nt] : died[i];
    for (const auto& interval : intervals) {
      for (int t = 0; t < nt; ++t) {
        const double lo = std::max({born[i], interval.first, time_breaks[t]});
        const double hi = std::min({life_end, interval.second, time_breaks[t + 1]});
        if (hi <= lo) continue;
        for (int a = 0; a < na; ++a) {
          const double enter = std::max(lo, born[i] + age_breaks[a]);
          const double exit = std::min(hi, born[i] + age_breaks[a + 1]);
          group_history.split(i, enter, exit, [&](int g, double lo, double hi) {
            partnered[t + nt * (a + na * s) + stride * g] += hi - lo;
          });
        }
      }
    }
  }
  NumericVector deaths_partnered(cells), deaths_unpartnered(cells);
  for (int i = 0; i < ids.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    const int t = demographic_bin(died[i], time_breaks);
    const int a = demographic_bin(died[i] - born[i], age_breaks);
    if (t < 0 || a < 0) continue;
    const int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
    const auto& intervals = episodes[i];
    // Left-limit status: the death closing an episode is partnered.
    const auto after = std::lower_bound(intervals.begin(), intervals.end(), died[i],
      [](const std::pair<double, double>& interval, double time) {
        return interval.first < time;
      });
    const bool in_union = after != intervals.begin() && died[i] <= std::prev(after)->second;
    if (in_union) deaths_partnered[t + nt * (a + na * s) + stride * group_history.at(i, died[i], true)] += 1;
    else deaths_unpartnered[t + nt * (a + na * s) + stride * group_history.at(i, died[i], true)] += 1;
  }
  NumericVector births_partnered(cells), births_unpartnered(cells);
  NumericVector mothers = persons["mother_id"];
  for (int child = 0; child < ids.size(); ++child) {
    if (child % 1024 == 0) checkUserInterrupt();
    const int t = demographic_bin(born[child], time_breaks);
    if (t < 0 || NumericVector::is_na(mothers[child])) continue;
    const int mother = index.at(mothers[child]);
    const int a = demographic_bin(born[child] - born[mother], age_breaks);
    if (a < 0) continue;
    const auto& intervals = episodes[mother];
    const auto after = std::upper_bound(intervals.begin(), intervals.end(), born[child],
      [](double time, const std::pair<double, double>& interval) {
        return time < interval.first;
      });
    const auto& dissolution_times = dissolution_ends[mother];
    const bool in_union = (after != intervals.begin() && born[child] < std::prev(after)->second) ||
      std::binary_search(dissolution_times.begin(), dissolution_times.end(), born[child]);
    if (in_union) births_partnered[t + nt * (a + na) + stride * group_history.at(mother, born[child], true)] += 1;
    else births_unpartnered[t + nt * (a + na) + stride * group_history.at(mother, born[child], true)] += 1;
  }
  return List::create(_["deaths_partnered"] = deaths_partnered,
                      _["deaths_unpartnered"] = deaths_unpartnered,
                      _["births_partnered"] = births_partnered,
                      _["births_unpartnered"] = births_unpartnered,
                      _["person_years_partnered"] = partnered);
}

// Union episodes, rather than persons: each record contributes exposure once.
// Death ends its exposure but only non-death dissolution enters the numerator.
// [[Rcpp::export]]
List dissolution_counts_cpp(DataFrame unions, NumericVector time_breaks,
                             NumericVector dur_breaks, IntegerVector group_indices,
                             int group_count, Nullable<NumericVector> origins = R_NilValue,
                             std::string event_reason = "dissolution") {
  GroupHistory group_history(unions, group_indices);
  NumericVector starts = unions["start_time"], ends = unions["end_time"];
  const NumericVector axis_origins = origins.isNotNull() ? NumericVector(origins) : starts;
  CharacterVector reasons = unions["end_reason"];
  const int nt = time_breaks.size() - 1, nd = dur_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * nd;
  NumericVector exposure(stride * group_count), events(stride * group_count);
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    const double end = NumericVector::is_na(ends[u]) ? time_breaks[nt] : ends[u];
    for (int t = 0; t < nt; ++t) {
      const double lo = std::max(starts[u], time_breaks[t]);
      const double hi = std::min(end, time_breaks[t + 1]);
      if (hi <= lo) continue;
      for (int d = 0; d < nd; ++d) {
        const double enter = std::max(lo, axis_origins[u] + dur_breaks[d]);
        const double exit = std::min(hi, axis_origins[u] + dur_breaks[d + 1]);
        group_history.split(u, enter, exit, [&](int g, double lo, double hi) {
          exposure[stride * g + t + nt * d] += hi - lo;
        });
      }
    }
    if (NumericVector::is_na(ends[u]) || CharacterVector::is_na(reasons[u]) ||
        as<std::string>(reasons[u]) != event_reason) continue;
    const int t = demographic_bin(ends[u], time_breaks);
    const int d = demographic_bin(ends[u] - axis_origins[u], dur_breaks);
    if (t >= 0 && d >= 0) events[stride * group_history.at(u, ends[u], true) + t + nt * d] += 1;
  }
  return List::create(_["union_years"] = exposure, _["dissolutions"] = events);
}

// Exclusive histories for previous-union state. Sort by time and union ID so
// input row order never changes status, including touching/zero-length episodes.
struct MarriageEpisode {
  double start, end, id;
  int after;
  bool cohab, conversion;
};
std::vector<std::vector<MarriageEpisode>> marriage_histories(DataFrame persons,
                                                            DataFrame unions) {
  NumericVector ids = persons["person_id"], died = persons["death_time"];
  NumericVector wives = unions["wife_id"], husbands = unions["husband_id"];
  NumericVector starts = unions["start_time"], ends = unions["end_time"], uids = unions["union_id"];
  CharacterVector reasons = unions["end_reason"];
  const bool typed = unions.containsElementNamed("union_type");
  CharacterVector types;
  if (typed) types = unions["union_type"];
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  std::vector<std::vector<MarriageEpisode>> histories(ids.size());
  for (int u = 0; u < starts.size(); ++u) {
    const bool ongoing = NumericVector::is_na(ends[u]);
    const double end = ongoing ? R_PosInf : ends[u];
    const int after = !ongoing && as<std::string>(reasons[u]) == "dissolution" ? 1 : 2;
    for (double id : {wives[u], husbands[u]}) {
      const int i = index.at(id);
      if (!NumericVector::is_na(died[i]) && (starts[u] > died[i] || end > died[i]))
        stop("Previous-union histories cannot contain union events after a partner's death.");
      histories[i].push_back({starts[u], end, uids[u], after,
        typed && as<std::string>(types[u]) == "cohab",
        !ongoing && as<std::string>(reasons[u]) == "conversion"});
    }
  }
  for (auto& history : histories) {
    std::sort(history.begin(), history.end(), [](const MarriageEpisode& a, const MarriageEpisode& b) {
      return a.start < b.start || (a.start == b.start && a.id < b.id);
    });
    for (size_t u = 1; u < history.size(); ++u)
      if (history[u].start < history[u - 1].end)
        stop("Previous-union stratification requires nonoverlapping union histories for each person.");
  }
  return histories;
}

// Execution-only cache: complete supplied history determines the latest ending.
// [[Rcpp::export]]
List previous_union_info_cpp(DataFrame persons, DataFrame unions, bool counts_cohab = false) {
  const auto histories = marriage_histories(persons, unions);
  IntegerVector state(histories.size(), 0), legal(histories.size(), 0);
  NumericVector ended(histories.size(), NA_REAL);
  for (size_t i = 0; i < histories.size(); ++i)
    for (const auto& episode : histories[i])
      if (R_finite(episode.end) && !episode.conversion) {
        if (!episode.cohab) legal[i] = episode.after;
        state[i] = episode.cohab && !counts_cohab ? legal[i] : episode.after;
        ended[i] = episode.end;
      }
  return List::create(_["previous_union"] = state, _["end_time"] = ended, _["legal_previous"] = legal);
}

// Existing internal snapshot interface for category-only validation.
// [[Rcpp::export]]
IntegerVector previous_union_state_cpp(DataFrame persons, DataFrame unions) {
  return as<IntegerVector>(previous_union_info_cpp(persons, unions)["previous_union"]);
}

// [[Rcpp::export]]
List marriage_counts_cpp(DataFrame persons, DataFrame unions,
                         NumericVector time_breaks, NumericVector age_breaks,
                         IntegerVector group_indices, int group_count,
                         Nullable<NumericVector> dur_breaks = R_NilValue, bool counts_cohab = false,
                         std::string entry_type = "pooled", bool marriage_measure = false) {
  GroupHistory group_history(persons, group_indices);
  const auto histories = marriage_histories(persons, unions);
  NumericVector born = persons["birth_time"], died = persons["death_time"];
  CharacterVector sex = persons["sex"];
  const int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const bool by_duration = dur_breaks.isNotNull();
  const NumericVector durations = by_duration ? NumericVector(dur_breaks) : age_breaks;
  const int nd = durations.size() - 1;
  const R_xlen_t age_cells = static_cast<R_xlen_t>(nt) * na * 2;
  const R_xlen_t dur_cells = static_cast<R_xlen_t>(nt) * nd * 2;
  const R_xlen_t stride = age_cells + 2 * dur_cells;
  const R_xlen_t offsets[] = {0, age_cells, age_cells + dur_cells};
  NumericVector exposure(stride * group_count), formations(exposure.size());
  for (int i = 0; i < born.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    const int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
    const double life_end = NumericVector::is_na(died[i]) ? time_breaks[nt] : died[i];
    auto add_exposure = [&](double lo, double hi, int previous) {
      const bool duration = by_duration && previous > 0;
      const NumericVector& breaks = duration ? durations : age_breaks;
      const int bands = breaks.size() - 1;
      const double origin = duration ? lo : born[i];
      for (int t = 0; t < nt; ++t) {
        const double start = std::max(lo, time_breaks[t]);
        const double end = std::min({hi, life_end, time_breaks[t + 1]});
        if (end <= start) continue;
        for (int a = 0; a < bands; ++a) {
          const double enter = std::max(start, origin + breaks[a]);
          const double exit = std::min(end, origin + breaks[a + 1]);
          group_history.split(i, enter, exit, [&](int g, double lo, double hi) {
            exposure[t + nt * (a + bands * s) + offsets[previous] + stride * g] += hi - lo;
          });
        }
      }
    };
    int previous = 0, legal = 0;
    bool converted_entry = false;
    double available = born[i];
    for (const auto& episode : histories[i]) {
      add_exposure(available, episode.start, previous);
      const int t = demographic_bin(episode.start, time_breaks);
      const bool duration = by_duration && previous > 0;
      const NumericVector& breaks = duration ? durations : age_breaks;
      const int bands = breaks.size() - 1;
      const double origin = duration ? available : born[i];
      const int a = demographic_bin(episode.start - origin, breaks);
      if ((marriage_measure ? !episode.cohab : !converted_entry) && t >= 0 && a >= 0 &&
          (entry_type == "pooled" || (episode.cohab ? entry_type == "cohab" : entry_type == "marriage")))
        formations[t + nt * (a + bands * s) + offsets[previous] + stride * group_history.at(i, episode.start, true)] += 1;
      if (marriage_measure && episode.cohab)
        add_exposure(episode.start, episode.end, previous);
      available = episode.end;
      converted_entry = episode.conversion;
      if (!episode.conversion) {
        if (!episode.cohab) legal = episode.after;
        previous = episode.cohab && !counts_cohab ? legal : episode.after;
      }
    }
    add_exposure(available, life_end, previous);
  }
  return List::create(_["person_years"] = exposure, _["unions"] = formations);
}

// Disjoint current-type demographic strata. Typed histories must be exclusive.
// [[Rcpp::export]]
List typed_demographic_counts_cpp(DataFrame persons, DataFrame unions,
    NumericVector time_breaks, NumericVector age_breaks,
    IntegerVector group_indices, int group_count, CharacterVector types) {
  GroupHistory group_history(persons, group_indices);
  const auto histories = marriage_histories(persons, unions);
  const int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * na * 2;
  const R_xlen_t cells = stride * group_count;
  NumericVector exposure(cells * types.size()), deaths(exposure.size()), births(exposure.size());
  CharacterVector union_types = unions["union_type"];
  Function subset("[.data.frame");
  for (int k = 0; k < types.size(); ++k) {
    LogicalVector rows(unions.nrows());
    for (int u = 0; u < unions.nrows(); ++u) rows[u] = union_types[u] == types[k];
    DataFrame selected = subset(unions, rows, R_MissingArg, Named("drop") = false);
    List x = union_counts_cpp(persons, selected, time_breaks, age_breaks, group_indices, group_count);
    NumericVector years = x["person_years_partnered"];
    std::copy(years.begin(), years.end(), exposure.begin() + cells * k);
  }
  NumericVector ids = persons["person_id"], born = persons["birth_time"], died = persons["death_time"], mothers = persons["mother_id"];
  CharacterVector sexes = persons["sex"];
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  auto classify = [&](int person, double time, bool birth) {
    const MarriageEpisode* active = nullptr;
    for (const auto& e : histories[person]) {
      // Birth precedes dissolution and conversion at exact ties.
      if (birth && e.end == time && (e.conversion || e.after == 1)) {
        active = &e;
        break;
      }
      if (birth ? (e.start <= time && time < e.end) : (e.start < time && time <= e.end))
        active = &e;
    }
    if (!active) return -1;
    const std::string label = active->cohab ? "cohab" : "marriage";
    for (int k = 0; k < types.size(); ++k) if (as<std::string>(types[k]) == label) return k;
    return -1;
  };
  for (int i = 0; i < ids.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    int t = demographic_bin(died[i], time_breaks), a = demographic_bin(died[i] - born[i], age_breaks);
    if (t >= 0 && a >= 0) {
      const int k = classify(i, died[i], false);
      const int sex = as<std::string>(sexes[i]) == "m" ? 0 : 1;
      if (k >= 0) deaths[cells * k + t + nt * (a + na * sex) + stride * group_history.at(i, died[i], true)] += 1;
    }
    t = demographic_bin(born[i], time_breaks);
    if (t < 0 || NumericVector::is_na(mothers[i])) continue;
    const int mother = index.at(mothers[i]);
    a = demographic_bin(born[i] - born[mother], age_breaks);
    if (a < 0) continue;
    const int k = classify(mother, born[i], true);
    if (k >= 0) births[cells * k + t + nt * (a + na) + stride * group_history.at(mother, born[i], true)] += 1;
  }
  return List::create(_["person_years"] = exposure, _["births"] = births, _["deaths"] = deaths);
}

// Row selection for parent-enriched father-choice observations. Merge union
// coverage so duplicate/overlapping episodes never duplicate births.
// [[Rcpp::export]]
IntegerVector unpartnered_birth_rows_cpp(DataFrame persons, DataFrame unions) {
  NumericVector ids = persons["person_id"], born = persons["birth_time"],
    mothers = persons["mother_id"], fathers = persons["father_id"];
  NumericVector wives = unions["wife_id"], starts = unions["start_time"], ends = unions["end_time"];
  CharacterVector reasons = unions["end_reason"];
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  std::vector<std::vector<std::pair<double, double>>> coverage(ids.size());
  std::vector<std::vector<double>> birth_first_ends(ids.size());
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    const int mother = index.at(wives[u]);
    coverage[mother].push_back({starts[u], NumericVector::is_na(ends[u]) ? R_PosInf : ends[u]});
    if (!NumericVector::is_na(ends[u]) && !CharacterVector::is_na(reasons[u])) {
      const std::string reason = as<std::string>(reasons[u]);
      if (reason == "dissolution" || reason == "conversion") birth_first_ends[mother].push_back(ends[u]);
    }
  }
  for (size_t i = 0; i < coverage.size(); ++i) {
    auto& intervals = coverage[i];
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<double, double>> merged;
    for (const auto& interval : intervals) {
      if (merged.empty() || interval.first > merged.back().second) merged.push_back(interval);
      else merged.back().second = std::max(merged.back().second, interval.second);
    }
    intervals = std::move(merged);
    std::sort(birth_first_ends[i].begin(), birth_first_ends[i].end());
  }
  std::vector<int> rows;
  for (int child = 0; child < ids.size(); ++child) {
    if (child % 1024 == 0) checkUserInterrupt();
    if (NumericVector::is_na(mothers[child]) || NumericVector::is_na(fathers[child])) continue;
    const int mother = index.at(mothers[child]);
    const auto& intervals = coverage[mother];
    const auto after = std::upper_bound(intervals.begin(), intervals.end(), born[child],
      [](double time, const std::pair<double, double>& interval) { return time < interval.first; });
    const auto& ties = birth_first_ends[mother];
    const bool partnered = (after != intervals.begin() && born[child] < std::prev(after)->second) ||
      std::binary_search(ties.begin(), ties.end(), born[child]);
    if (!partnered) rows.push_back(child + 1);
  }
  return wrap(rows);
}

// Observed maternal fertility: split lifetimes at births before Lexis clipping.
// [[Rcpp::export]]
List fertility_parity_counts_cpp(DataFrame persons, DataFrame unions,
    NumericVector time_breaks, NumericVector age_breaks, NumericVector parity_breaks,
    IntegerVector group_indices, int group_count, CharacterVector types, bool split_types) {
  GroupHistory group_history(persons, group_indices);
  NumericVector ids=persons["person_id"], born=persons["birth_time"],
    died=persons["death_time"], mothers=persons["mother_id"];
  CharacterVector sexes=persons["sex"];
  int nt=time_breaks.size()-1, na=age_breaks.size()-1, np=parity_breaks.size()-1;
  int nk=types.size();
  R_xlen_t base=static_cast<R_xlen_t>(nt)*na*group_count;
  NumericVector exposure(base*np*nk), events(exposure.size());
  std::map<double,int> index;
  for(int i=0;i<ids.size();++i) index[ids[i]]=i;
  std::vector<std::vector<double>> births(ids.size());
  for(int i=0;i<ids.size();++i) if(!NumericVector::is_na(mothers[i]))
    births[index.at(mothers[i])].push_back(born[i]);
  struct Episode { double start,end; int type; bool birth_first; };
  std::vector<std::vector<Episode>> histories(ids.size());
  if(split_types) {
    NumericVector wife=unions["wife_id"], starts=unions["start_time"], ends=unions["end_time"];
    CharacterVector reasons=unions["end_reason"];
    bool typed=unions.containsElementNamed("union_type");
    CharacterVector labels=typed ? as<CharacterVector>(unions["union_type"]) : CharacterVector(unions.nrows(),"partnered");
    for(int u=0;u<unions.nrows();++u) {
      int k=1;
      for(;k<nk;++k) if(types[k]==labels[u]) break;
      bool first=!CharacterVector::is_na(reasons[u]) &&
        (reasons[u]=="dissolution" || reasons[u]=="conversion");
      histories[index.at(wife[u])].push_back({starts[u],NumericVector::is_na(ends[u])?R_PosInf:ends[u],k,first});
    }
  }
  for(int i=0;i<ids.size();++i) {
    if(i%1024==0) checkUserInterrupt();
    if(sexes[i]!="f") continue;
    auto &b=births[i]; std::stable_sort(b.begin(),b.end());
    auto &h=histories[i];
    std::sort(h.begin(),h.end(),[](const Episode&a,const Episode&b){return a.start<b.start;});
    // Merge same-type coverage for legacy overlapping episodes.
    std::vector<Episode> coverage;
    for(auto e:h) {
      if(!coverage.empty() && e.type==coverage.back().type && e.start<=coverage.back().end)
        coverage.back().end=std::max(coverage.back().end,e.end);
      else coverage.push_back(e);
    }
    auto cell=[&](int t,int a,int p,int k,int g){return t+nt*a+static_cast<R_xlen_t>(nt)*na*g+base*(p+np*k);};
    double life_end=NumericVector::is_na(died[i])?time_breaks[nt]:died[i];
    for(size_t j=0;j<=b.size();++j) {
      int p=demographic_bin(static_cast<double>(j),parity_breaks);
      double left=j?b[j-1]:born[i], right=j<b.size()?b[j]:life_end;
      if(p<0) continue;
      for(int t=0;t<nt;++t) for(int a=0;a<na;++a) {
        double lo=std::max({left,born[i],time_breaks[t],born[i]+age_breaks[a]});
        double hi=std::min({right,life_end,time_breaks[t+1],born[i]+age_breaks[a+1]});
        if(hi<=lo) continue;
        group_history.split(i, lo, hi, [&](int g, double lo, double hi) {
          double remaining=hi-lo;
          if(split_types) for(const auto&e:coverage) {
            double years=std::max(0.0,std::min(hi,e.end)-std::max(lo,e.start));
            exposure[cell(t,a,p,e.type,g)]+=years; remaining-=years;
          }
          exposure[cell(t,a,p,0,g)]+=std::max(0.0,remaining);
        });
      }
      if(j<b.size()) {
        int t=demographic_bin(b[j],time_breaks), a=demographic_bin(b[j]-born[i],age_breaks), k=0;
        if(split_types) {
          for(const auto&e:h) if(e.start<=b[j] && b[j]<e.end) k=e.type;
          for(const auto&e:h) if(e.end==b[j] && e.birth_first) { k=e.type; break; }
        }
        if(t>=0 && a>=0) events[cell(t,a,p,k,group_history.at(i,b[j],true))]+=1;
      }
    }
  }
  return List::create(_["person_years"]=exposure,_["births"]=events);
}

// Directed transitions: origin exposure is stored once, not per destination.
// [[Rcpp::export]]
List transition_counts_cpp(DataFrame persons, DataFrame unions, DataFrame transitions,
    IntegerVector final_groups, int group_count, NumericVector time_breaks,
    NumericVector axis_breaks, bool duration, bool by_sex, bool by_union_type) {
  NumericVector born = persons["birth_time"], died = persons["death_time"];
  CharacterVector sexes = persons["sex"];
  const int nt = time_breaks.size() - 1, na = axis_breaks.size() - 1;
  const int ns = by_sex ? 2 : 1;
  const bool typed = unions.containsElementNamed("union_type");
  const int nk = !by_union_type ? 1 : (typed ? 3 : 2);
  const R_xlen_t base = static_cast<R_xlen_t>(nt) * na;
  NumericVector exposure(base * group_count * ns * nk);
  NumericVector events(base * group_count * group_count * ns * nk);
  struct Change { double time; int from, to; };
  std::vector<std::vector<Change>> changes(persons.nrows());
  IntegerVector row = transitions["row"], from = transitions["from"], to = transitions["to"];
  NumericVector times = transitions["time"];
  for (int j = 0; j < times.size(); ++j) changes[row[j]].push_back({times[j], from[j], to[j]});
  std::vector<std::vector<MarriageEpisode>> unions_by_person;
  if (by_union_type) unions_by_person = marriage_histories(persons, unions);
  for (int i = 0; i < persons.nrows(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    auto& moves = changes[i];
    std::sort(moves.begin(), moves.end(), [](const Change& a, const Change& b) { return a.time < b.time; });
    const int sex = by_sex && sexes[i] == "f" ? 1 : 0;
    auto union_type = [&](double time) {
      if (!by_union_type) return 0;
      const auto& episodes = unions_by_person[i];
      auto next = std::upper_bound(episodes.begin(), episodes.end(), time,
        [](double t, const MarriageEpisode& e) { return t < e.start; });
      if (next == episodes.begin()) return 0;
      const auto& active = *std::prev(next);
      if (time >= active.end) return 0;
      return typed ? (active.cohab ? 1 : 2) : 1;
    };
    auto exposure_cell = [&](int t, int a, int g, int k) -> R_xlen_t {
      return t + static_cast<R_xlen_t>(nt) * a + base * (g + static_cast<R_xlen_t>(group_count) * (sex + ns * k));
    };
    auto event_cell = [&](int t, int a, int origin, int dest, int k) -> R_xlen_t {
      return t + static_cast<R_xlen_t>(nt) * a + base * (dest + static_cast<R_xlen_t>(group_count) *
        (origin + static_cast<R_xlen_t>(group_count) * (sex + ns * k)));
    };
    int origin = moves.empty() ? final_groups[i] : moves.front().from;
    double entry = born[i];
    const double life_end = NumericVector::is_na(died[i]) ? time_breaks[nt] : died[i];
    for (size_t j = 0; j <= moves.size(); ++j) {
      const double end = j < moves.size() ? moves[j].time : life_end;
      const double lo = std::max(entry, time_breaks[0]);
      const double hi = std::min({end, life_end, time_breaks[nt]});
      const double axis_origin = duration ? entry : born[i];
      if (hi > lo) {
        std::vector<double> cuts = {lo, hi};
        if (by_union_type) for (const auto& e : unions_by_person[i]) {
          if (lo < e.start && e.start < hi) cuts.push_back(e.start);
          if (lo < e.end && e.end < hi) cuts.push_back(e.end);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        for (size_t c = 1; c < cuts.size(); ++c) {
          const int k = union_type(cuts[c - 1]);
          for (int t = 0; t < nt; ++t) {
            const double start = std::max(cuts[c - 1], time_breaks[t]);
            const double stop = std::min(cuts[c], time_breaks[t + 1]);
            if (stop <= start) continue;
            for (int a = 0; a < na; ++a) {
              const double enter = std::max(start, axis_origin + axis_breaks[a]);
              const double exit = std::min(stop, axis_origin + axis_breaks[a + 1]);
              if (exit > enter) exposure[exposure_cell(t, a, origin, k)] += exit - enter;
            }
          }
        }
      }
      if (j < moves.size()) {
        const auto& move = moves[j];
        const int t = demographic_bin(move.time, time_breaks);
        const int a = demographic_bin(move.time - axis_origin, axis_breaks);
        if (t >= 0 && a >= 0) events[event_cell(t, a, move.from, move.to, union_type(move.time))] += 1;
        origin = move.to;
        entry = move.time;
      }
    }
  }
  return List::create(_["person_years"] = exposure, _["transitions"] = events);
}

// Partner identity at formation, before a same-time group transition. Reuse
// indexed histories to avoid repeated R scans and named-list lookups per union.
// [[Rcpp::export]]
IntegerMatrix union_start_groups_cpp(DataFrame persons, DataFrame unions,
                                     IntegerVector group_indices) {
  GroupHistory history(persons, group_indices);
  NumericVector ids = persons["person_id"], wives = unions["wife_id"],
    husbands = unions["husband_id"], starts = unions["start_time"];
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  IntegerMatrix groups(unions.nrows(), 2);
  for (int u = 0; u < unions.nrows(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    groups(u, 0) = history.at(index.at(wives[u]), starts[u], true);
    groups(u, 1) = history.at(index.at(husbands[u]), starts[u], true);
  }
  return groups;
}
