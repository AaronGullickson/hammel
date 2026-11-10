# Kinship module and generation-distance exclusions

Specified October 15, 2026 from the user's request for shared native kinship
code, a two-generation default, and screening only demonstrated relationships.
**Status:** implemented recorded-parent kinship for union screening and R extraction.
No genetic relatedness coefficient, social-kin inference, or legal-marital rule.

## Scientific rule and worked examples

Let D(x,a) be the minimum number of recorded parent steps from person x to
ancestor a. Include the person themself at distance zero. At configured limit G,
exclude a candidate man y for woman x if any recorded common ancestor a
satisfies BOTH `D(x,a) <= G` and `D(y,a) <= G`. This is a pair of generation
bounds, not a bound on their sum or an undirected family-graph distance.

The user selected `kinship_generations = 2` as the default. It is a finite,
nonnegative whole number. Zero imposes no kin exclusion between distinct
people; two includes known first cousins. Increasing the limit also increases
the range of direct ancestry examined. There is no separate unlimited exclusion
of direct ancestors. Minimum male age remains independently configurable.

| Ego distance | Relative distance | Simple-pedigree example | Excluded at G = 2? |
| --- | --- | --- | --- |
| 1 | 0 | Parent of ego | Yes |
| 0 | 1 | Child of ego | Yes |
| 1 | 1 | Sibling, including a half sibling | Yes |
| 2 | 0 | Grandparent of ego | Yes |
| 2 | 1 | Aunt/uncle of ego | Yes |
| 1 | 2 | Niece/nephew of ego | Yes |
| 2 | 2 | First cousin, including a half first cousin | Yes |
| 2 | 3 | First cousin once removed | No; included at G = 3 |
| 3 | 3 | Second cousin | No; included at G = 3 |
| 3 | 0 | Great-grandparent of ego | No; included at G = 3 |

Only demonstrated links count. Unknown mothers/fathers are separate missing
links, never a shared artificial ancestor. A known paternal relationship is
still screened when a mother is unknown. No detected relationship does not
prove unrelatedness. Unrelated co-parents do not become blood kin through a
shared child. Death does not remove an ancestor or its recorded links.

The exact partner draw (D072) checks minimum male age and kinship for each
drawn man, redrawing when he fails. If no man in the pool passes, birth attempts are
cancelled by default (D043); `cancel_birth_if_no_match = FALSE` retains a birth
with unknown father and no new union. Failed rate searches reschedule. Supplied existing unions are retained;
these are restrictions on new searches, not retrospective history corrections.
Kinship adds no minimum female age. Independently, D037 requires exclusive
ongoing unions; partnered men are not in the available pool.

## Shared native implementation contract

Use a separate `kinship.h`/`kinship.cpp` module, called directly from the engine
and through an Rcpp extraction entry point. Resolve public IDs to internal
positions at the boundary. Parent links are immutable; append a newborn's
recorded parents before scheduling its future events. Rebuild the same graph
on restart. All graph/cache state belongs to the call.

Validate references and reject ancestry cycles before querying. Do not infer
parent sex/age rules beyond existing population validation. Iterate rather than
recursively consuming the C++ stack. Keep deceased persons as ancestry nodes.

An ancestry query performs breadth-first traversal up to G, visiting each
DISTINCT ancestor once at its minimum distance. Cache queried people's bounded
ancestor sets. For screening, mark the current woman's ancestor membership
once and check candidates by direct lookup, stopping on the first match. Clear
only the previous woman's marked entries when the query ego changes; retain
sorted sets for R extraction. This is one call-local membership vector, not
a dense all-pairs matrix. At a zero bound, skip candidate kinship queries and
ancestry-cache building entirely; input cycle validation remains.
Retain maternal/paternal/both lineage information for shortest routes. Reusable
visit marks avoid allocating a hash entry for every route. Cache validity for
existing people is unaffected by newborn append because their parents do not
change. No dense all-pairs matrix or enumeration of every ancestral path.

At the default depth of two, a cached set has at most seven entries (self,
two parents, and four grandparents). Traversal can still grow when ancestors
are distinct at deeper limits. Work
and stored cache entries are bounded by reachable distinct people, not by
materializing repeated routes; depth and memory costs remain explicit. Full
path enumeration and genetic coefficients are outside this increment.

## R extraction interface

`extract_kin(population, ego_ids, generations = 2, relative_ids = NULL)` returns
an ordinary tibble of known relatives within the SAME two-sided bounds. Egos
are explicitly selected; relative IDs default to all people in recorded history.
Accept a `hammel_population` object with a finite reference time (D166).
Bare persons tables must first be constructed into a population with an explicit
`end_time`. Both deceased and living relatives are eligible for extraction. This is a
history query, not a living-kin census at a selected date (see
`count_kin()` below). Inputs must be validated populations with
valid, unique query IDs.

It returns ONE row per ego/relative pair (D161). Columns: `ego_id`,
`relative_id`, `relationship`, `ego_distance`, `relative_distance`,
`n_common_ancestors`, `ego_lineage`, `relative_lineage`,
`ego_ancestry_complete`, and `relative_ancestry_complete`. Completeness means
all parent links needed WITHIN the requested depth are recorded; it says
nothing about more distant generations. Missing ancestry does not suppress
known evidence. No rows means no demonstrated relationship within bounds,
not proof of unrelatedness or complete ancestry.

The closest common ancestors are those with the smallest SUM of the two
distances among the bounded shared ancestors. The row reports their distance
pair and `n_common_ancestors`, the number of closest recorded ancestors at that
pair: ordinarily 2 for full kin (both parents, or two grandparents) and 1 for
half kin. A count of 1 is only recorded evidence: it does not prove half kin
when the other parent is unknown. Lineages combine across those ancestors:
`maternal` and `paternal` routes together give `both`. In rare collapsed
pedigrees, equally close ancestors can sit at different distance pairs; the
row keeps the pair with the smaller ego distance, and counts only the
ancestors at that pair. This output convention is for compact extraction
only: screening rejects ANY shared ancestor within both bounds.

`relationship` names the distance pair (ego distance, relative distance), with
full and half kin together:

| Pair | Label | Pair | Label |
| --- | --- | --- | --- |
| (1,0) | `parent` | (0,1) | `child` |
| (2,0) | `grandparent` | (0,2) | `grandchild` |
| (k,0), k >= 3 | `great_grandparent`, `great_great_grandparent`, ... | (0,k) | `great_grandchild`, ... |
| (1,1) | `sibling` | (2,2) | `first_cousin` |
| (2,1) | `aunt_uncle` | (1,2) | `niece_nephew` |
| (k,1), k >= 3 | `great_aunt_uncle`, ... | (1,k) | `great_niece_nephew`, ... |
| (a,b), a, b >= 2 | `<n>_cousin`, n = min(a,b) - 1 | removal | `_once_removed`, `_twice_removed`, ... for |a - b| = 1, 2, ... |

For example (3,3) is `second_cousin` and (2,3) and (3,2) are both
`first_cousin_once_removed`. Ordinals run `first` to `fifth`, then `6th`;
removals run `once`, `twice`, `thrice`, then `4_times`.

Use an ancestor-to-selected-relatives index to find extraction candidates
instead of comparing every ego with every person. Preserve input ego and
relative order. Empty results have the same typed tibble columns. No candidate
search sets are retained.

## Living kin counts

`count_kin(population, times, generations = 2, by = NULL)` averages
living kin per living person at each requested time (D161). Inputs follow
`calculate_population_size()`: a `hammel_population` whose known `end_time` is
no earlier than the requested finite times. A person is living at time t when
`birth_time <= t` and `death_time` is missing or after t, as for population
size: births at t count and deaths at t do not.

At each time, every living person is an ego. A relative counts when the two
share a recorded ancestor within `generations` on each side (as in
`extract_kin()`) and the relative is also living at t. Each relative counts
once, under the relationship of the closest shared ancestors; full and half
kin are pooled. Kin outside the bound are ignored.

`by` takes `"sex"` and declared group dimensions, describing the EGO at t
(after events at t, as in `calculate_population_size()`). Output is a tibble
with `time`, the `by` columns, `relationship`, `persons` (living egos in the
stratum), `kin` (their living kin of that type, summed), `mean_kin`
(`kin / persons`, NA when there are no egos) and `share_with_kin` (share of
egos with at least one, NA when there are no egos). Every time, stratum and
relationship possible within `generations` has a row, including zero counts.
Distance pairs that share a label, such as (2,3) and (3,2) for
`first_cousin_once_removed`, form one relationship. Rows vary by time fastest,
then the `by` variables in order, then relationship in increasing distance sum (parent, child, grandparent,
sibling, grandchild, aunt_uncle, niece_nephew, first_cousin at depth two).

No standard errors are reported: relatives share kin, so egos are not
independent, and replicate simulations are the appropriate source of
uncertainty. Ancestry completeness is not used to filter egos. Counts reflect
recorded links only, so a starting population without recorded parents
undercounts kin until two or three generations have been simulated; users
must allow a burn-in before interpreting counts.

Worked example (also a test): grandparents 1 (f) and 2 (m) born -60, with 1
dying at -10; their children 3 (f) and 4 (m) born -30; unrelated 5 (f) and
6 (m) born -30. Person 7 (born -5) has parents 3 and 6; person 8 (born -5)
has parents 5 and 4; person 9 (born -5) has mother 3 and an unknown father.
For ego 7, `extract_kin()` gives grandparents 1 and 2 (2,0) with one ancestor
each, parents 3 and 6 (1,0), aunt/uncle 4 (2,1) with two ancestors, first
cousin 8 (2,2) with two ancestors and half sibling 9 (1,1) with one. At time
0, ego 7 has living kin: 1 grandparent (2), 2 parents, 1 sibling, 1
aunt/uncle and 1 first cousin. At time -20 only persons 1 to 6 are living
(1 dies later, at -10): 3 and 4 each have 2 living parents and 1 sibling, 1
and 2 each have 2 children, and 5 and 6 have no recorded kin. Over those 6
egos, `mean_kin` is 4/6 for parents, 4/6 for children and 2/6 for siblings,
and 0 for the other types.

## Performance

`count_kin()` extracts each person's bounded kin once, then counts
living egos and kin at each time with vector operations. Time grows linearly
with the number of recorded persons and about threefold from depth two to
three. On November 9, 2026 (commit `fc60749` plus these changes; optimized
install, R 4.6.1, 10-core ARM Linux), 11 times over 100 simulated years from
USA 2000 rates took (median of three runs; `by = "sex"` adds nothing
measurable):

| Starting persons | Recorded persons | Depth 2 | Depth 3 |
| --- | --- | --- | --- |
| 5,000 | 12,175 | 0.14 s | 0.44 s |
| 20,000 | 48,261 | 0.64 s | 1.72 s |
| 50,000 | 122,210 | 1.87 s | 5.39 s |
| 100,000 | 244,231 | 4.53 s | 12.7 s |

Group dimensions with recorded transitions add a per-person history lookup at
each time, which this benchmark does not exercise. Rerun with
`validation/kinship/benchmark-count-kin.R`.

## Validation

Use hand-built maternal/paternal pedigrees with noncontiguous and shuffled IDs.
Cover all worked examples, half kin, missing parents, unrelated co-parents,
multiple routes/pedigree collapse, deceased ancestors, empty schemas, invalid
IDs/depths, cycles, and the depth-zero limit. Compare extraction and engine
screening on the same pairs. Verify dynamically appended children cannot
be selected by their own mothers, continued-run kin evidence matches rebuilt
histories, and independent calls do not share state. Scientific model recovery
and genetic/social kin extensions are separate work.
