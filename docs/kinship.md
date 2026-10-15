# Kinship module and generation-distance exclusions

Specified October 15, 2026 from the user's request for shared native kinship
code, a two-generation default, and screening only demonstrated relationships.
**Status:** implemented recorded-parent kinship for union screening and R extraction.
No genetic relatedness coefficient, social-kin inference, or legal-marital rule.

## Scientific rule and worked examples

Let D(x,a) be the minimum number of recorded parent steps from person x to
ancestor a. Include the person themself at distance zero. At configured limit G,
exclude a sampled candidate y for woman x if any recorded common ancestor a
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

Screen minimum male age and kinship AFTER sampling and BEFORE matching-model
scoring. Do not refill the sample. If no candidate passes, retain the birth
with unknown father and no new union. Supplied existing unions are retained;
these are restrictions on new searches, not retrospective history corrections.
No extra minimum female age or male exclusivity is introduced.

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
Accept a population object or a base persons data frame. No known endpoint
is required. Both deceased and living relatives are eligible for extraction. This is a
history query, not an automatic living-kin census at a selected date. Inputs
must be validated populations with valid, unique query IDs.

Columns: `ego_id`, `relative_id`, `common_ancestor_id`, `ego_distance`,
`relative_distance`, `ego_lineage`, `relative_lineage`,
`ego_ancestry_complete`, and `relative_ancestry_complete`. Completeness means
all parent links needed WITHIN the requested depth are recorded; it says
nothing about more distant generations. Missing ancestry does not suppress
known evidence. No rows means no demonstrated relationship within bounds,
not proof of unrelatedness or complete ancestry.

For each ego/relative pair, return common ancestors achieving the smallest
SUM of the two distances among the bounded shared ancestors. Return ties
rather than invent a unique lineage. This output convention is for compact
extraction only: screening rejects ANY shared ancestor within both bounds.
Lineage is `self`, `maternal`, `paternal`, or `both`, describing the first
parental side(s) on shortest routes to the reported ancestor. Do not infer
full versus half relationships from absent evidence.

Thus full siblings ordinarily produce two rows at (1,1), one for each shared
parent, rather than also being listed as (2,2) through their grandparents.
First cousins can have two equally close common grandparents. More complex
histories can have multiple equally close routes/relationships. No categorical
relationship label is assigned; distances/evidence avoid silently resolving
ambiguous categories. Count distinct relative IDs when counting people.

Use an ancestor-to-selected-relatives index to find extraction candidates
instead of comparing every ego with every person. Preserve input ego and
relative order, sorting tied ancestor IDs for stable output. Empty results
have the same typed tibble columns. No candidate search sets are retained.

## Validation

Use hand-built maternal/paternal pedigrees with noncontiguous and shuffled IDs.
Cover all worked examples, half kin, missing parents, unrelated co-parents,
multiple routes/pedigree collapse, deceased ancestors, empty schemas, invalid
IDs/depths, cycles, and the depth-zero limit. Compare extraction and engine
screening on the same pairs. Verify dynamically appended children cannot
be selected by their own mothers, continued-run kin evidence matches rebuilt
histories, and independent calls do not share state. Scientific model recovery
and genetic/social kin extensions are separate work.
