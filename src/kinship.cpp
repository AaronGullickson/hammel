#include <Rcpp.h>
#include "kinship.h"
#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace hammel {

Kinship::Kinship(std::vector<std::array<int, 2>> parents, int generations)
  : generations_(generations), parents_(std::move(parents)), cache_(parents_.size()),
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
        complete = false;
        continue;
      }
      const int distance = current.distance + 1;
      const unsigned char lineage = current.distance == 0 ?
        static_cast<unsigned char>(1 << side) : current.lineage;
      if (visited_[parent] != stamp_) {
        visited_[parent] = stamp_;
        positions_[parent] = static_cast<int>(result.size());
        result.push_back({parent, distance, lineage});
      } else {
        auto& prior = result[positions_[parent]];
        if (prior.distance == distance) prior.lineage |= lineage;
      }
    }
  }
  std::sort(result.begin(), result.end(), [](const KinAncestor& a, const KinAncestor& b) {
    return a.person < b.person;
  });
  ready_[person] = 1;
  complete_[person] = complete;
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
  std::vector<double> out_ego, out_relative, out_ancestor;
  std::vector<int> out_ego_distance, out_relative_distance;
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
      auto evidence = kinship.shared_ancestors(ego, relative);
      long long nearest = std::numeric_limits<long long>::max();
      for (const auto& shared : evidence)
        nearest = std::min(nearest, static_cast<long long>(shared.ego_distance) + shared.relative_distance);
      std::sort(evidence.begin(), evidence.end(), [&](const hammel::KinEvidence& a, const hammel::KinEvidence& b) {
        return ids[a.ancestor] < ids[b.ancestor];
      });
      for (const auto& shared : evidence) {
        if (static_cast<long long>(shared.ego_distance) + shared.relative_distance != nearest) continue;
        out_ego.push_back(ego_ids[e]); out_relative.push_back(relative_ids[r]);
        out_ancestor.push_back(ids[shared.ancestor]);
        out_ego_distance.push_back(shared.ego_distance);
        out_relative_distance.push_back(shared.relative_distance);
        out_ego_lineage.push_back(kin_lineage(shared.ego_lineage));
        out_relative_lineage.push_back(kin_lineage(shared.relative_lineage));
        out_ego_complete.push_back(kinship.ancestry_complete(ego));
        out_relative_complete.push_back(kinship.ancestry_complete(relative));
      }
    }
  }
  return DataFrame::create(_["ego_id"] = out_ego, _["relative_id"] = out_relative,
    _["common_ancestor_id"] = out_ancestor, _["ego_distance"] = out_ego_distance,
    _["relative_distance"] = out_relative_distance, _["ego_lineage"] = out_ego_lineage,
    _["relative_lineage"] = out_relative_lineage,
    _["ego_ancestry_complete"] = LogicalVector(out_ego_complete.begin(), out_ego_complete.end()),
    _["relative_ancestry_complete"] = LogicalVector(out_relative_complete.begin(), out_relative_complete.end()));
}
