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
  const int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  const R_xlen_t stride = static_cast<R_xlen_t>(nt) * na * 2;
  const R_xlen_t cells = stride * group_count;
  NumericVector partnered(cells), formations(cells);
  std::map<double, int> index;
  std::vector<std::vector<std::pair<double, double>>> episodes(ids.size());
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  for (int u = 0; u < starts.size(); ++u) {
    if (u % 1024 == 0) checkUserInterrupt();
    for (double id : {wives[u], husbands[u]}) {
      const int i = index.at(id);
      const int s = as<std::string>(sex[i]) == "m" ? 0 : 1;
      const double end = NumericVector::is_na(ends[u]) ? time_breaks[nt] : ends[u];
      episodes[i].push_back({starts[u], end});
      const int t = demographic_bin(starts[u], time_breaks);
      const int a = demographic_bin(starts[u] - born[i], age_breaks);
      if (t >= 0 && a >= 0) formations[t + nt * (a + na * s) + stride * group_indices[i]] += 1;
    }
  }
  for (int i = 0; i < ids.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
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
    const bool in_union = after != intervals.begin() && born[child] < std::prev(after)->second;
    if (in_union) births_partnered[t + nt * (a + na) + stride * group_indices[mother]] += 1;
    else births_unpartnered[t + nt * (a + na) + stride * group_indices[mother]] += 1;
  }
  return List::create(_["deaths_partnered"] = deaths_partnered,
                      _["deaths_unpartnered"] = deaths_unpartnered,
                      _["births_partnered"] = births_partnered,
                      _["births_unpartnered"] = births_unpartnered,
                      _["person_years_partnered"] = partnered,
                      _["unions"] = formations);
}
