#ifndef HAMMEL_KINSHIP_H
#define HAMMEL_KINSHIP_H

#include <array>
#include <cstdint>
#include <vector>

namespace hammel {

// Positions in the caller's persons vector, not public person IDs. -1 is unknown.
struct KinAncestor {
  int person;
  int distance;
  unsigned char lineage; // 0 self, 1 maternal, 2 paternal, 3 both shortest sides
};

struct KinEvidence {
  int ancestor;
  int ego_distance;
  int relative_distance;
  unsigned char ego_lineage;
  unsigned char relative_lineage;
};

class Kinship {
public:
  Kinship(std::vector<std::array<int, 2>> parents, int generations);
  void append(int mother, int father);
  int generations() const { return generations_; }
  const std::vector<KinAncestor>& ancestors(int person);
  bool known_related(int first, int second);
  std::vector<KinEvidence> shared_ancestors(int ego, int relative);
  bool ancestry_complete(int person);

private:
  int generations_;
  std::vector<std::array<int, 2>> parents_;
  std::vector<std::vector<KinAncestor>> cache_;
  std::vector<unsigned char> ready_, complete_;
  std::vector<std::uint64_t> visited_;
  std::vector<int> positions_;
  std::uint64_t stamp_ = 0;
  // A single ego's ancestor membership, not a population-by-population matrix.
  std::vector<unsigned char> query_marks_;
  std::vector<int> marked_ancestors_;
  int query_person_ = -1;
  void validate_acyclic() const;
  void check_person(int person) const;
};

} // namespace hammel
#endif
