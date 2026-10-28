#include <Rcpp.h>
#include <algorithm>
#include <map>

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
        if (exit > enter) exposure[t + nt * (a + na * s) + stride * group_indices[i]] += exit - enter;
      }
    }
    int t = demographic_bin(died[i], time_breaks);
    int a = demographic_bin(died[i] - born[i], age_breaks);
    if (t >= 0 && a >= 0) deaths[t + nt * (a + na * s) + stride * group_indices[i]] += 1;

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
    if (a >= 0) births[t + nt * (a + na) + stride * group_indices[mother->second]] += 1;
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
          if (exit > enter) partnered[t + nt * (a + na * s) + stride * group_indices[i]] += exit - enter;
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
    if (in_union) deaths_partnered[t + nt * (a + na * s) + stride * group_indices[i]] += 1;
    else deaths_unpartnered[t + nt * (a + na * s) + stride * group_indices[i]] += 1;
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
    if (in_union) births_partnered[t + nt * (a + na) + stride * group_indices[mother]] += 1;
    else births_unpartnered[t + nt * (a + na) + stride * group_indices[mother]] += 1;
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
                             int group_count) {
  NumericVector starts = unions["start_time"], ends = unions["end_time"];
  CharacterVector reasons = unions["end_reason"];
  const int nt = time_breaks.size() - 1, nd = dur_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * nd;
  NumericVector exposure(stride * group_count), events(stride * group_count);
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    const R_xlen_t offset = stride * group_indices[u];
    const double end = NumericVector::is_na(ends[u]) ? time_breaks[nt] : ends[u];
    for (int t = 0; t < nt; ++t) {
      const double lo = std::max(starts[u], time_breaks[t]);
      const double hi = std::min(end, time_breaks[t + 1]);
      if (hi <= lo) continue;
      for (int d = 0; d < nd; ++d) {
        const double enter = std::max(lo, starts[u] + dur_breaks[d]);
        const double exit = std::min(hi, starts[u] + dur_breaks[d + 1]);
        if (exit > enter) exposure[offset + t + nt * d] += exit - enter;
      }
    }
    if (NumericVector::is_na(ends[u]) || CharacterVector::is_na(reasons[u]) ||
        as<std::string>(reasons[u]) != "dissolution") continue;
    const int t = demographic_bin(ends[u], time_breaks);
    const int d = demographic_bin(ends[u] - starts[u], dur_breaks);
    if (t >= 0 && d >= 0) events[offset + t + nt * d] += 1;
  }
  return List::create(_["union_years"] = exposure, _["dissolutions"] = events);
}

// Exclusive histories for previous-union state. Sort by time and union ID so
// input row order never changes status, including touching/zero-length episodes.
struct MarriageEpisode {
  double start, end, id;
  int after;
};
std::vector<std::vector<MarriageEpisode>> marriage_histories(DataFrame persons,
                                                            DataFrame unions) {
  NumericVector ids = persons["person_id"], died = persons["death_time"];
  NumericVector wives = unions["wife_id"], husbands = unions["husband_id"];
  NumericVector starts = unions["start_time"], ends = unions["end_time"], uids = unions["union_id"];
  CharacterVector reasons = unions["end_reason"];
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
      histories[i].push_back({starts[u], end, uids[u], after});
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
IntegerVector previous_union_state_cpp(DataFrame persons, DataFrame unions) {
  const auto histories = marriage_histories(persons, unions);
  IntegerVector state(histories.size(), 0);
  for (size_t i = 0; i < histories.size(); ++i)
    for (const auto& episode : histories[i])
      if (R_finite(episode.end)) state[i] = episode.after;
  return state;
}

// [[Rcpp::export]]
List marriage_counts_cpp(DataFrame persons, DataFrame unions,
                         NumericVector time_breaks, NumericVector age_breaks,
                         IntegerVector group_indices, int group_count) {
  const auto histories = marriage_histories(persons, unions);
  NumericVector born = persons["birth_time"], died = persons["death_time"];
  CharacterVector sex = persons["sex"];
  const int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * na * 2;
  NumericVector exposure(stride * 3 * group_count), formations(exposure.size());
  for (int i = 0; i < born.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    const int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
    const double life_end = NumericVector::is_na(died[i]) ? time_breaks[nt] : died[i];
    auto add_exposure = [&](double lo, double hi, int previous) {
      for (int t = 0; t < nt; ++t) {
        const double start = std::max(lo, time_breaks[t]);
        const double end = std::min({hi, life_end, time_breaks[t + 1]});
        if (end <= start) continue;
        for (int a = 0; a < na; ++a) {
          const double enter = std::max(start, born[i] + age_breaks[a]);
          const double exit = std::min(end, born[i] + age_breaks[a + 1]);
          if (exit > enter)
            exposure[t + nt * (a + na * s) + stride * (previous + 3 * group_indices[i])] += exit - enter;
        }
      }
    };
    int previous = 0;
    double available = born[i];
    for (const auto& episode : histories[i]) {
      add_exposure(available, episode.start, previous);
      const int t = demographic_bin(episode.start, time_breaks);
      const int a = demographic_bin(episode.start - born[i], age_breaks);
      if (t >= 0 && a >= 0)
        formations[t + nt * (a + na * s) + stride * (previous + 3 * group_indices[i])] += 1;
      available = episode.end;
      previous = episode.after;
    }
    add_exposure(available, life_end, previous);
  }
  return List::create(_["person_years"] = exposure, _["unions"] = formations);
}
