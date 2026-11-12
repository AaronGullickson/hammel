# External immigration: design direction

Planning discussion, November 12, 2026. External immigration is not
implemented. D177–D182 record the accepted direction; the open questions below
must be resolved before implementation. Illustrative field names do not
establish a public API. See the [roadmap](roadmap.md#post-release-immigration-from-outside-the-system)
for scope. No emigration mechanism is selected here.

## Arrival process and streams

Use an absolute arrival intensity in migrants per year (D177), rather than a
hazard applied to residents. For individual arrivals at constant intensity
lambda, exponential interarrival times produce a Poisson count over T years,
with mean and variance lambda * T. An intensity of 100 implies an expected
100 migrants in one year and a standard deviation of 10, not a guaranteed
count. Resident population changes do not themselves change the intensity.
Controlled/exact-count immigration is excluded.

Allow multiple streams (D178), potentially representing countries or world
regions of origin, each with its own intensity and source template. Stream
identity does not automatically become current group membership, ancestry or
an inherited origin category. Calendar-varying intensity inputs and clock
organization still need specification.

## Population templates and selection

Draw migrant templates with replacement from an existing `hammel_population`
(D179). Preserve age, sex and group identifiers, and assign fresh destination
person IDs on every draw. Preserve age by translating birth time from the
source population's endpoint to the arrival date. The template is a reusable
snapshot and is not depleted by migration.

The user proposed a predictive selection model as an alternative to helpers
that tune the source population's distributions (D180). The proposed
separation is that arrival intensity determines the expected flow, while a
selection model determines who is sampled conditional on an arrival.
Normalized relative weights would preserve the absolute stream intensity:
multiplying every eligible person's weight by the same positive constant
would change neither the flow nor the selection distribution. This remains a
proposed interpretation. Relative weights, predicted probabilities and
individual hazards are not interchangeable. Source eligibility, model inputs
and outputs, age effects, zero weights and empty eligible pools still need
decisions. The existing `weight_model()` contract cannot simply be assumed
to cover this new use.

Allow each stream's template to change over calendar time using time bounds
and a list column of population objects (D180). Column naming, validation,
coverage and boundary semantics remain open. Selecting the active template
when an arrival occurs, with changes affecting future arrivals only, is the
proposed event-time rule.

Optional household immigration is under discussion. Reuse of household
transition membership rules, household selection probabilities, within-party
ID/link reconstruction and relatives outside the party need specification.
The user explicitly deferred whether intensity means households per year or
migrants per year. D177's individual-arrival example does not settle the
person-count distribution for household arrivals.

## Migration facts, family history and classifications

Keep migration history separate from derived generation/origin classifications
(D181). One generation number cannot describe every mixed-stream pedigree.

| Layer | Meaning |
| --- | --- |
| Personal migration facts | Whether a person actually immigrated into the destination, when, and through which stream |
| Family migration history | Connections to migration histories through parent–child relationships |
| Derived classifications | Named variables under a specified rule, usable in demographic mechanisms and analysis |

For example, a locally born child may have a mother who immigrated through
stream A and a locally born father whose parents immigrated through stream B.
The child has an immigrant parent connected to A and immigrant grandparents
connected to B. Preserve enough information to distinguish those connections;
do not silently select one origin or a universal generation rule.

An immigrant arriving as a child has their own arrival history. A source
relative retained only to support an arriving person's reproductive history
has not thereby immigrated. Source genealogy alone does not establish a
relative's migration history with respect to the destination.

Migration facts, group membership and the package's ancestry proportions are
distinct. Classification may use several streams without requiring a group
category for each stream. Mappings among these concepts must be explicit.

## Concise access in custom demographic models

Users should be able to define a migration classification once and reuse its
named outputs in partner matching, group inheritance, weight models and
analysis (D181). Custom functions should receive derived values in their
documented inputs rather than having to traverse pedigrees. This is a new
interface requirement, not existing callback functionality.

Illustrative outputs include an immigrant indicator, personal arrival stream,
a generation summary under a chosen rule, and origin-specific generation
classifications. These are examples, not selected field names, default
generation rules or a mandatory fixed set of variables. A simple model may
need one summary, while a mixed-origin model may need several outputs.

Specify classification inputs, output names and types, unknown-history
behavior, and initialization/update timing, especially at birth and
immigration. Distinguish unknown history from a known absence of connection
to a stream. Custom mechanisms must retain existing state-consistency and
event-invalidation requirements; unrestricted state-changing callbacks are
not authorized.

## Rate stratification

Support rate schedules stratified by personal immigration status and arrival
stream (D182), without requiring a generational classification in the rate
table. Illustrative categories distinguish people who have not themselves
immigrated, immigrants through stream A, and immigrants through stream B.
These are conceptual categories, not selected column names or encodings.

A destination-born child of stream-A immigrants belongs to the category of
people who have not themselves immigrated. Their family migration history
still connects them to stream A for D181 classifications and custom models.
Arrival-stream status records personal immigration, not inherited origin.
Unknown starting-population histories must not silently be classified as a
known absence of immigration; their handling remains to be specified.

Status/stream-specific schedules set aggregate rate targets for those strata.
Derived generational variables in normalized weight models instead influence
who experiences events within rate cells, without changing cell intensities.
They do not independently impose generation-specific aggregate rates.

Whether to use an ordinary declared group or a built-in migration stratifier
remains open, as do process coverage, schedule validation and category names.
General generation-specific rate stratification is not selected. This design
does not change the current implementation's requirement that categorical
rate strata use declared group dimensions.

## Truncation and history retention

The requested truncation helper retains living people and their descendants,
including deceased descendants (D180). This preserves recorded births needed
for parity and last-birth reconstruction; it does not impute unknown births.
Deceased ancestors need not be retained solely because living people descend
from them. See the [burn-in helper discussion](roadmap.md#living-population-extraction-after-burn-in).

Removing ancestors can discard migration information needed for classification.
Preserve sufficient history or summaries before removal so supported
classifications do not change simply because records were pruned. The
representation and supported history depth remain unresolved. Preserving
today's classification label alone may be insufficient for future descendants
or for applying a different rule.

Truncation also needs rules for removed parent references, current/previous
unions, ancestry, observation windows and restarting. Retaining selected
historical records does not imply complete burn-in histories or pedigrees.
Immigration separately needs arrival-aware exposure and event accounting so
imported pre-arrival births and histories are not mistaken for destination
events or exposure.

## Specification and validation still needed

Settle history/attribute copying, group and ancestry compatibility, generation
rules and representation, selection semantics, schedules, arrival/observation
timing, event ties, household behavior and public constructors before
implementation. No implementation milestone or scientific acceptance threshold
is established here.

Validation should follow those choices. Relevant cases include replicated
arrival counts and waiting times, independent-stream scenarios, replacement
with fresh IDs, age preservation, schedule boundaries, mixed-stream pedigrees,
classification inputs in custom models, and classification/parity preservation
through truncation and restart. Include unknown histories and distinguish
actual arrivals from supporting historical relatives. Household checks depend
on its deferred sampling and intensity decisions.
