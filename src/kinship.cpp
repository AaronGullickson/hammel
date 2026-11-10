#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace hammel {

Kinship::Kinship(std::vector<std::array<int, 2>> parents, int generations,
                 bool track_details)
  : generations_(generations), track_details_(track_details),
    parents_(std::move(parents)), cache_(parents_.size()),
    ready_(parents_.size(), 0), complete_(parents_.size(), 0),
    visited_(parents_.size(), 0), positions_(parents_.size(), 0),
    query_marks_(parents_.size(), 0) {
  if (generations < 0) throw std::invalid_argument("Kinship generations must be nonnegative.");
  for (const auto& links : parents_) {
    for (int parent : links) {
      if (parent < -1 || (parent >= 0 && static_cast<size_t>(parent) >= parents_.size()))
        throw std::invalid_argument("Kinship parent reference is invalid.");
    }
  }
  validate_acyclic();
}

void Kinship::check_person(int person) const {
  if (person < 0 || static_cast<size_t>(person) >= parents_.size())
    throw std::invalid_argument("Kinship person reference is invalid.");
}

void Kinship::validate_acyclic() const {
  // Iterative DFS: 0 unseen, 1 on current ancestor path, 2 fully checked.
  std::vector<unsigned char> state(parents_.size(), 0);
  std::vector<std::pair<int, int>> stack;
  for (int person = 0; person < static_cast<int>(parents_.size()); ++person) {
    if (state[person]) continue;
    state[person] = 1;
    stack.emplace_back(person, 0);
    while (!stack.empty()) {
      auto& frame = stack.back();
      if (frame.second == 2) {
        state[frame.first] = 2;
        stack.pop_back();
        continue;
      }
      const int parent = parents_[frame.first][frame.second++];
      if (parent < 0) continue;
      if (state[parent] == 1)
        throw std::invalid_argument("Kinship parent links contain an ancestry cycle.");
      if (state[parent] == 0) {
        state[parent] = 1;
        stack.emplace_back(parent, 0);
      }
    }
  }
}

void Kinship::append(int mother, int father) {
  // Newborn parents already exist, so appending cannot introduce an ancestry cycle.
  for (int parent : {mother, father}) {
    if (parent != -1) check_person(parent);
  }
  parents_.push_back({mother, father});
  cache_.emplace_back();
  ready_.push_back(0);
  complete_.push_back(0);
  visited_.push_back(0);
  positions_.push_back(0);
  query_marks_.push_back(0);
}

const std::vector<KinAncestor>& Kinship::ancestors(int person) {
  check_person(person);
  if (ready_[person]) return cache_[person];
  if (++stamp_ == 0) {
    std::fill(visited_.begin(), visited_.end(), 0);
    stamp_ = 1;
  }
  auto& result = cache_[person];
  result.push_back({person, 0, 0});
  visited_[person] = stamp_;
  positions_[person] = 0;
  bool complete = true;
  // Breadth-first traversal retains minimum distances and only distinct ancestors.
  for (size_t cursor = 0; cursor < result.size(); ++cursor) {
    const KinAncestor current = result[cursor];
    if (current.distance >= generations_) continue;
    for (int side = 0; side < 2; ++side) {
      const int parent = parents_[current.person][side];
      if (parent < 0) {
        if (track_details_) {
          complete = false;
        }
        continue;
      }
      const int distance = current.distance + 1;
      const unsigned char lineage = !track_details_ ? 0 :
        (current.distance == 0 ?
          static_cast<unsigned char>(1 << side) : current.lineage);
      if (visited_[parent] != stamp_) {
        visited_[parent] = stamp_;
        positions_[parent] = static_cast<int>(result.size());
        result.push_back({parent, distance, lineage});
      } else if (track_details_) {
        auto& prior = result[positions_[parent]];
        if (prior.distance == distance) prior.lineage |= lineage;
      }
    }
  }
  std::sort(result.begin(), result.end(), [](const KinAncestor& a, const KinAncestor& b) {
    return a.person < b.person;
  });
  ready_[person] = 1;
  if (track_details_) {
    complete_[person] = complete;
  }
  return result;
}

bool Kinship::known_related(int first, int second) {
  // The same woman is checked against a batch of candidates. Mark her ancestry
  // once, then scan only each candidate's cached set. Clear only the previous
  // marked set when the ego changes, without clearing a population-sized array.
  if (query_person_ != first) {
    const auto& ego_ancestors = ancestors(first);
    for (int person : marked_ancestors_) query_marks_[person] = 0;
    marked_ancestors_.clear();
    for (const auto& ancestor : ego_ancestors) {
      query_marks_[ancestor.person] = 1;
      marked_ancestors_.push_back(ancestor.person);
    }
    query_person_ = first;
  }
  for (const auto& ancestor : ancestors(second)) {
    if (query_marks_[ancestor.person]) return true;
  }
  return false;
}

std::vector<KinEvidence> Kinship::shared_ancestors(int ego, int relative) {
  if (!track_details_) {
    throw std::logic_error("Kinship reporting details were not requested.");
  }
  const auto& a = ancestors(ego);
  const auto& b = ancestors(relative);
  std::vector<KinEvidence> result;
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i].person == b[j].person) {
      result.push_back({a[i].person, a[i].distance, b[j].distance,
                        a[i].lineage, b[j].lineage});
      ++i; ++j;
    } else if (a[i].person < b[j].person) ++i;
    else ++j;
  }
  return result;
}

bool Kinship::ancestry_complete(int person) {
  if (!track_details_) {
    throw std::logic_error("Kinship reporting details were not requested.");
  }
  ancestors(person);
  return complete_[person];
}

} // namespace hammel

using namespace Rcpp;

static const char* kin_lineage(unsigned char lineage) {
  switch (lineage) {
  case 1: return "maternal";
  case 2: return "paternal";
  case 3: return "both";
  default: return "self";
  }
}

// [[Rcpp::export]]
DataFrame extract_kin_cpp(DataFrame persons, NumericVector ego_ids,
                          NumericVector relative_ids, int generations) {
  NumericVector ids = persons["person_id"], mothers = persons["mother_id"];
  NumericVector fathers = persons["father_id"];
  std::map<double, int> index;
  for (int i = 0; i < ids.size(); ++i) index[ids[i]] = i;
  std::vector<std::array<int, 2>> parents(ids.size());
  for (int i = 0; i < ids.size(); ++i) {
    parents[i] = {NumericVector::is_na(mothers[i]) ? -1 : index.at(mothers[i]),
                  NumericVector::is_na(fathers[i]) ? -1 : index.at(fathers[i])};
  }
  hammel::Kinship kinship(std::move(parents), generations);
  // Invert the bounded ancestor sets for selected relatives, not all person pairs.
  std::vector<std::vector<int>> relatives_of_ancestor(ids.size());
  for (int r = 0; r < relative_ids.size(); ++r) {
    if (r % 1024 == 0) checkUserInterrupt();
    const int relative = index.at(relative_ids[r]);
    for (const auto& ancestor : kinship.ancestors(relative))
      relatives_of_ancestor[ancestor.person].push_back(r);
  }
  std::vector<double> out_ego, out_relative;
  std::vector<int> out_ego_distance, out_relative_distance, out_ancestors;
  std::vector<std::string> out_ego_lineage, out_relative_lineage;
  std::vector<int> out_ego_complete, out_relative_complete;
  std::vector<int> seen(relative_ids.size(), -1);
  for (int e = 0; e < ego_ids.size(); ++e) {
    checkUserInterrupt();
    const int ego = index.at(ego_ids[e]);
    std::vector<int> candidates;
    for (const auto& ancestor : kinship.ancestors(ego)) {
      for (int r : relatives_of_ancestor[ancestor.person]) {
        if (seen[r] != e) {
          seen[r] = e;
          if (relative_ids[r] != ego_ids[e]) candidates.push_back(r);
        }
      }
    }
    std::sort(candidates.begin(), candidates.end());
    for (int r : candidates) {
      const int relative = index.at(relative_ids[r]);
      const auto evidence = kinship.shared_ancestors(ego, relative);
      // Closest ancestors minimize the distance sum; among them, keep the
      // pair with the smaller ego distance (differs only in collapsed
      // pedigrees). One row per pair counts and combines those ancestors.
      int best_sum = std::numeric_limits<int>::max();
      int best_ego = std::numeric_limits<int>::max();
      for (const auto& shared : evidence) {
        const int sum = shared.ego_distance + shared.relative_distance;
        if (sum < best_sum || (sum == best_sum && shared.ego_distance < best_ego)) {
          best_sum = sum;
          best_ego = shared.ego_distance;
        }
      }
      int count = 0;
      unsigned char ego_lineage = 0, relative_lineage = 0;
      for (const auto& shared : evidence) {
        if (shared.ego_distance != best_ego ||
            shared.ego_distance + shared.relative_distance != best_sum) continue;
        ++count;
        ego_lineage |= shared.ego_lineage;
        relative_lineage |= shared.relative_lineage;
      }
      if (!count) continue;
      out_ego.push_back(ego_ids[e]); out_relative.push_back(relative_ids[r]);
      out_ego_distance.push_back(best_ego);
      out_relative_distance.push_back(best_sum - best_ego);
      out_ancestors.push_back(count);
      out_ego_lineage.push_back(kin_lineage(ego_lineage));
      out_relative_lineage.push_back(kin_lineage(relative_lineage));
      out_ego_complete.push_back(kinship.ancestry_complete(ego));
      out_relative_complete.push_back(kinship.ancestry_complete(relative));
    }
  }
  return DataFrame::create(_["ego_id"] = out_ego, _["relative_id"] = out_relative,
    _["ego_distance"] = out_ego_distance,
    _["relative_distance"] = out_relative_distance,
    _["n_common_ancestors"] = out_ancestors, _["ego_lineage"] = out_ego_lineage,
    _["relative_lineage"] = out_relative_lineage,
    _["ego_ancestry_complete"] = LogicalVector(out_ego_complete.begin(), out_ego_complete.end()),
    _["relative_ancestry_complete"] = LogicalVector(out_relative_complete.begin(), out_relative_complete.end()));
}

// Count selected relationships without allocating a table of extracted pairs.
// `strata` is person x observation time: zero means not living, otherwise a
// one-based ego stratum. Relationship lookup entries are zero-based types,
// or -1 for excluded distance pairs. Output varies time, stratum, then type.
// [[Rcpp::export]]
List count_kin_cpp(DataFrame persons, int generations,
                   IntegerMatrix relationship_types, IntegerMatrix strata,
                   int n_strata, int n_types) {
  NumericVector ids = persons["person_id"], mothers = persons["mother_id"];
  NumericVector fathers = persons["father_id"];
  const int n = ids.size();
  const int n_times = strata.ncol();
  if (generations < 0 || strata.nrow() != n || n_strata < 0 || n_types < 0 ||
      relationship_types.nrow() != generations + 1 ||
      relationship_types.ncol() != generations + 1) {
    stop("Invalid native kin-count dimensions.");
  }
  for (int type : relationship_types) {
    if (type < -1 || type >= n_types) {
      stop("Invalid native kin-count relationship index.");
    }
  }
  std::map<double, int> index;
  for (int i = 0; i < n; ++i) {
    index[ids[i]] = i;
  }
  std::vector<std::array<int, 2>> parents(n);
  for (int i = 0; i < n; ++i) {
    parents[i] = {
      NumericVector::is_na(mothers[i]) ? -1 : index.at(mothers[i]),
      NumericVector::is_na(fathers[i]) ? -1 : index.at(fathers[i])
    };
  }
  // Preserve whole-graph cycle validation, even for zero requested types.
  hammel::Kinship kinship(std::move(parents), generations, false);
  const R_xlen_t n_cells = static_cast<R_xlen_t>(n_times) * n_strata;
  IntegerVector persons_n(n_cells), kin_n(n_cells * n_types);
  IntegerVector with_n(n_cells * n_types);
  std::vector<int> active;
  for (int person = 0; person < n; ++person) {
    bool living = false;
    for (int t = 0; t < n_times; ++t) {
      const int stratum = strata(person, t);
      if (stratum < 0 || stratum > n_strata) {
        stop("Invalid native kin-count stratum.");
      }
      if (stratum == 0) {
        continue;
      }
      living = true;
      ++persons_n[t + static_cast<R_xlen_t>(n_times) * (stratum - 1)];
    }
    if (living) {
      active.push_back(person);
    }
  }
  if (n_types == 0) {
    return List::create(_["persons"] = persons_n, _["kin"] = kin_n,
                        _["with"] = with_n);
  }
  struct Relative {
    int person;
    int distance;
  };
  // Retain only ancestor postings with a requested relative-side distance.
  std::vector<unsigned char> relative_distances(generations + 1, 0);
  for (int b = 0; b <= generations; ++b) {
    for (int a = 0; a <= generations; ++a) {
      if (relationship_types(a, b) >= 0) {
        relative_distances[b] = 1;
      }
    }
  }
  std::vector<std::vector<Relative>> relatives_of_ancestor(n);
  for (int relative : active) {
    if (relative % 1024 == 0) {
      checkUserInterrupt();
    }
    for (const auto& ancestor : kinship.ancestors(relative)) {
      if (relative_distances[ancestor.distance]) {
        relatives_of_ancestor[ancestor.person].push_back(
          {relative, ancestor.distance});
      }
    }
  }
  std::vector<int> seen(n, -1);
  std::vector<unsigned char> has_kin(
    static_cast<R_xlen_t>(n_times) * n_types, 0);
  for (int ego : active) {
    checkUserInterrupt();
    const auto& ego_ancestors = kinship.ancestors(ego);
    std::vector<int> candidates;
    for (const auto& ancestor : ego_ancestors) {
      for (const auto& relative : relatives_of_ancestor[ancestor.person]) {
        if (relative.person == ego || seen[relative.person] == ego ||
            relationship_types(ancestor.distance, relative.distance) < 0) {
          continue;
        }
        seen[relative.person] = ego;
        candidates.push_back(relative.person);
      }
    }
    std::fill(has_kin.begin(), has_kin.end(), 0);
    for (int relative : candidates) {
      const auto& other = kinship.ancestors(relative);
      int best_sum = std::numeric_limits<int>::max();
      int best_ego = std::numeric_limits<int>::max();
      size_t a = 0, b = 0;
      // Always examine the full ancestry bound, including unselected closer
      // relationships. Selection must not turn siblings into cousins.
      while (a < ego_ancestors.size() && b < other.size()) {
        if (ego_ancestors[a].person == other[b].person) {
          const int distance = ego_ancestors[a].distance;
          const int sum = distance + other[b].distance;
          if (sum < best_sum || (sum == best_sum && distance < best_ego)) {
            best_sum = sum;
            best_ego = distance;
          }
          ++a;
          ++b;
        } else if (ego_ancestors[a].person < other[b].person) {
          ++a;
        } else {
          ++b;
        }
      }
      const int type = relationship_types(best_ego, best_sum - best_ego);
      if (type < 0) {
        continue;
      }
      for (int t = 0; t < n_times; ++t) {
        const int stratum = strata(ego, t);
        if (stratum == 0 || strata(relative, t) == 0) {
          continue;
        }
        const R_xlen_t cell = t + static_cast<R_xlen_t>(n_times) *
          (stratum - 1 + static_cast<R_xlen_t>(n_strata) * type);
        if (kin_n[cell] == std::numeric_limits<int>::max()) {
          stop("Kin counts exceed integer range.");
        }
        ++kin_n[cell];
        has_kin[t + static_cast<R_xlen_t>(n_times) * type] = 1;
      }
    }
    for (int type = 0; type < n_types; ++type) {
      for (int t = 0; t < n_times; ++t) {
        if (has_kin[t + static_cast<R_xlen_t>(n_times) * type]) {
          const R_xlen_t cell = t + static_cast<R_xlen_t>(n_times) *
            (strata(ego, t) - 1 + static_cast<R_xlen_t>(n_strata) * type);
          ++with_n[cell];
        }
      }
    }
  }
  return List::create(_["persons"] = persons_n, _["kin"] = kin_n,
                      _["with"] = with_n);
}
