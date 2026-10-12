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
                             NumericVector age_breaks) {
  NumericVector ids = persons["person_id"], born = persons["birth_time"];
  NumericVector died = persons["death_time"], mothers = persons["mother_id"];
  CharacterVector sex = persons["sex"];
  int nt = time_breaks.size() - 1, na = age_breaks.size() - 1;
  R_xlen_t cells = static_cast<R_xlen_t>(nt) * na * 2;
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
        if (exit > enter) exposure[t + nt * (a + na * s)] += exit - enter;
      }
    }
    int t = demographic_bin(died[i], time_breaks);
    int a = demographic_bin(died[i] - born[i], age_breaks);
    if (t >= 0 && a >= 0) deaths[t + nt * (a + na * s)] += 1;

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
    if (a >= 0) births[t + nt * (a + na)] += 1;
  }
  return List::create(_["person_years"] = exposure, _["deaths"] = deaths,
                       _["births"] = births, _["unassigned_births"] = unknown);
}
