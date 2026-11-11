// Piecewise-constant rate schedules (bands) read from R rate tables, and
// selection of the calendar period active at a time.
#ifndef HAMMEL_RATE_BANDS_H
#define HAMMEL_RATE_BANDS_H

#include <Rcpp.h>

#include <string>
#include <utility>
#include <vector>

namespace hammel {

using namespace Rcpp;

// One piecewise constant age interval: [start, end), in years, with an annual
// person-year rate. An open final interval has end = infinity.
struct Band {
  double start, end, rate;
};

// Copy the rows of an R rate table that match a selection into bands sorted
// by start, without changing the R input. Empty strings and negative values
// select all rows: `sex` ("m" or "f"), `status` (union type), `group` (the
// .hammel_group index, when present), `previous` (previous-union category)
// and `parity` (rows whose parity band holds it). `duration` reads duration
// bands instead of age bands.
std::vector<Band> read_bands(
  DataFrame rates,
  const std::string& sex = "",
  const std::string& status = "",
  int group = -1,
  const std::string& previous = "",
  bool duration = false,
  double parity = -1
);

// Select active flat-table rows at a calendar boundary. The public table is
// never changed. R validation checks every active slice before execution.
DataFrame calendar_slice(DataFrame rates, double now);

// Whether any calendar period of a rate table starts or ends at now.
bool calendar_changes_at(DataFrame rates, double now);

// Whether search for a previous-union category runs on a duration clock
// (time since the last union ended) rather than age. R validation guarantees
// one clock per previous category across all groups.
bool search_uses_duration(DataFrame rates, const std::string& previous);

// Actual current hazard and next absolute age boundary, including schedule
// gaps. Absolute comparisons make a queued boundary use the right-hand band
// even when subtracting birth time would round the attained age just below its
// edge.
inline std::pair<double, double>
current_hazard(const std::vector<Band>* bands, double born, double now) {
  double value = 0, next = R_PosInf;
  if (bands) {
    for (const Band& band : *bands) {
      const double start = born + band.start, end = born + band.end;
      if (start <= now && now < end) {
        value = band.rate;
      }
      if (start > now) {
        next = std::min(next, start);
      }
      if (end > now) {
        next = std::min(next, end);
      }
    }
  }
  return {value, next};
}

// Locate time `now` on sorted, disjoint bands measured from `origin`: the
// index of the band holding it (-1 in a gap or outside every band), with the
// next absolute band edge in `next` (infinity if none). Absolute comparisons,
// as in current_hazard().
inline int locate_band(
  const std::vector<Band>& bands,
  double origin,
  double now,
  double& next
) {
  size_t lo = 0, hi = bands.size();
  while (lo < hi) {
    const size_t mid = (lo + hi) / 2;
    if (origin + bands[mid].start <= now) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  next = R_PosInf;
  if (lo > 0 && now < origin + bands[lo - 1].end) {
    next = origin + bands[lo - 1].end;
    return static_cast<int>(lo - 1);
  }
  if (lo < bands.size()) {
    next = origin + bands[lo].start;
  }
  return -1;
}

}  // namespace hammel

#endif  // HAMMEL_RATE_BANDS_H
