// Fertility helpers declared in fertility_cells.h.
#include "fertility_cells.h"

#include <Rcpp.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "rate_bands.h"

using namespace Rcpp;

namespace hammel {

double
inherit_fertility_weight(double maternal, double variance, double alpha) {
  if (alpha == 1) {
    return maternal;
  }
  const double mu = -variance / 2;
  double log_weight = mu + alpha * (std::log(maternal) - mu);
  if (variance > 0) {
    log_weight +=
      std::sqrt(variance) * std::sqrt(1 - alpha * alpha) * R::rnorm(0, 1);
  }
  const double weight = std::exp(log_weight);
  if (!(weight > 0) || !R_finite(weight)) {
    stop("Generated hetfert weight is not finite and strictly positive.");
  }
  return weight;
}

std::vector<FertilityCells::Source> fertility_strata(
  DataFrame rates,
  int groups,
  const std::vector<std::string>& type_labels
) {
  std::vector<FertilityCells::Source> out;
  const bool parity = rates.containsElementNamed("parity_start");
  const bool grouped = rates.containsElementNamed(".hammel_group");
  for (int g = 0; g < groups; ++g) {
    for (int t = 0; t < static_cast<int>(type_labels.size()); ++t) {
      const std::string label = type_labels[t];
      std::map<double, double> bands = {{0, R_PosInf}};
      if (parity) {
        bands.clear();
        NumericVector starts = rates["parity_start"],
                      ends = rates["parity_end"];
        IntegerVector row_groups;
        if (grouped) {
          row_groups = rates[".hammel_group"];
        }
        CharacterVector row_types;
        if (!label.empty()) {
          row_types = rates["union_type"];
        }
        for (int i = 0; i < starts.size(); ++i) {
          if (grouped && row_groups[i] != g) {
            continue;
          }
          if (!label.empty() && as<std::string>(row_types[i]) != label) {
            continue;
          }
          bands[starts[i]] = ends[i];
        }
      }
      for (const auto& band : bands) {
        FertilityCells::Source st;
        st.group = g;
        st.type = t;
        st.parity_start = band.first;
        st.parity_end = band.second;
        st.bands =
          read_bands(rates, "", label, g, "", false, parity ? band.first : -1);
        out.push_back(st);
      }
    }
  }
  return out;
}

}  // namespace hammel
