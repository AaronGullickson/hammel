// Rate-table readers declared in rate_bands.h.
#include "rate_bands.h"

#include <Rcpp.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace Rcpp;

namespace hammel {

std::vector<Band> read_bands(
  DataFrame rates,
  const std::string& sex,
  const std::string& status,
  int group,
  const std::string& previous,
  bool duration,
  double parity
) {
  NumericVector starts = rates[duration ? "dur_start" : "age_start"],
                ends = rates[duration ? "dur_end" : "age_end"];
  NumericVector values = rates["rate"];
  CharacterVector sexes;
  if (!sex.empty()) {
    sexes = rates["sex"];
  }
  CharacterVector statuses;
  if (!status.empty()) {
    statuses = rates["union_type"];
  }
  CharacterVector previous_types;
  if (!previous.empty()) {
    previous_types = rates["previous_union"];
  }
  IntegerVector groups;
  const bool grouped = rates.containsElementNamed(".hammel_group");
  if (grouped) {
    groups = rates[".hammel_group"];
  }
  NumericVector parity_start, parity_end;
  const bool select_parity =
    parity >= 0 && rates.containsElementNamed("parity_start");
  if (select_parity) {
    parity_start = rates["parity_start"];
    parity_end = rates["parity_end"];
  }
  std::vector<Band> bands;
  for (int i = 0; i < starts.size(); ++i) {
    if (select_parity &&
        !(parity_start[i] <= parity && parity < parity_end[i])) {
      continue;
    }
    if ((sex.empty() || as<std::string>(sexes[i]) == sex) &&
        (status.empty() || as<std::string>(statuses[i]) == status) &&
        (!grouped || groups[i] == group) &&
        (previous.empty() || as<std::string>(previous_types[i]) == previous)) {
      bands.push_back({starts[i], ends[i], values[i]});
    }
  }
  std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) {
    return a.start < b.start;
  });
  return bands;
}

DataFrame calendar_slice(DataFrame rates, double now) {
  if (!rates.containsElementNamed("time_start")) {
    return rates;
  }
  NumericVector lo = rates["time_start"], hi = rates["time_end"];
  std::vector<int> rows;
  for (int i = 0; i < rates.nrows(); ++i) {
    if (lo[i] <= now && now < hi[i]) {
      rows.push_back(i + 1);
    }
  }
  Function subset = Environment::base_env()["["];
  return as<DataFrame>(
    subset(rates, wrap(rows), R_MissingArg, _["drop"] = false)
  );
}

bool calendar_changes_at(DataFrame rates, double now) {
  if (!rates.containsElementNamed("time_start")) {
    return false;
  }
  NumericVector lo = rates["time_start"], hi = rates["time_end"];
  for (int i = 0; i < rates.nrows(); ++i) {
    if (lo[i] == now || hi[i] == now) {
      return true;
    }
  }
  return false;
}

bool search_uses_duration(DataFrame rates, const std::string& previous) {
  if (!rates.containsElementNamed("dur_start") || previous.empty()) {
    return false;
  }
  CharacterVector labels = rates["previous_union"];
  NumericVector ages = rates["age_start"];
  for (int i = 0; i < labels.size(); ++i) {
    if (as<std::string>(labels[i]) == previous) {
      return NumericVector::is_na(ages[i]);
    }
  }
  return false;
}

}  // namespace hammel
