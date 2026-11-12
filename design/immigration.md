# External immigration design

Planning discussion, November 12, 2026. External immigration is not
implemented. This document consolidates the approach selected in D177–D182
for design review on the `external-migration` branch. Behavior depending on
the unresolved modeling choices below requires further user decisions before
implementation. Illustrative field names do not establish a public API.
See the [roadmap](roadmap.md#post-release-immigration-from-outside-the-system)
for scope. No emigration mechanism is selected here.

## Purpose and intended workflow

Add people from outside the simulated population during a continuous-time
run. Users should be able to describe different origin populations and
migration flows without constructing every migrant record manually or
restarting the destination simulation for each arrival.

The intended workflow is:

1. Prepare a source `hammel_population`, potentially from a separate burn-in
   simulation. A template need not be stable or represent the whole origin
   population; users control the supplied source.
2. Optionally truncate source histories while preserving living people and
   their descendants, including deceased descendants needed for reproductive
   histories.
3. Define one or more streams, each referencing a template and specifying an
   absolute arrival intensity. Templates may change through calendar schedules.
4. Optionally describe selective migration from the template. A predictive
   selection model is the current proposal; its exact contract is open.
5. Run the destination model, record actual arrivals and preserve sufficient
   migration history for rate strata, classifications and later analysis.

Population preparation, stream/model construction and simulation execution
remain separate responsibilities. No constructor names or signatures are
selected by this workflow.

## Decision summary

| Topic | Current position | Record |
| --- | --- | --- |
| Arrival flow | Absolute stochastic intensity; no controlled/exact-count option | D177 |
| Origins | Multiple streams with separate intensities | D178 |
| Source | Sample an existing population with replacement; preserve age, sex and groups; assign fresh IDs | D179 |
| Template preparation and schedules | Retain living people and descendants when truncating; support calendar-scheduled templates | D180 |
| Selective migration | Predictive selection proposed instead of distribution-tuning helpers; model semantics open | D180 |
| Migration histories | Separate personal/family facts from reusable, named classifications | D181 |
| Rate schedules | Stratify by personal immigration status and stream; generation rules stay outside the required rate-table schema | D182 |
| Households | Optional household migration under discussion; sampling and intensity units deferred | D179–D180 |

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

## Integration with the destination simulation

Arrivals must compete chronologically with the existing continuous-time
processes. An immigration event adds new destination members rather than
selecting an existing resident to experience an event. Whether streams have
separate clocks or share a clock with stream selection is not decided.

An arriving person must enter all applicable destination process cells and
matching pools with consistent relationships and demographic state. Source
event clocks cannot be transplanted: destination rates and models determine
subsequent events. The initialization of parity, spacing, unions, group
duration and other history-dependent state still needs a scientific contract.
Importing supporting histories is distinct from treating every retained source
relative as an arriving member of the destination population.

Migration records and outputs should remain ordinary, inspectable R objects,
with public tables following the package's tibble convention. Their schema,
ID mapping/provenance and restart representation remain open. Each run must
own its state and must not mutate a reusable source template as it samples.

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

## Open decisions

The following choices are open; listing an option does not select it.

| Area | Decisions still needed |
| --- | --- |
| Source sampling | Eligible source records; uniform versus modeled selection; weight/probability interpretation; empty templates and zero eligible weight |
| Imported state | Which histories and attributes accompany a draw; relatives outside a party; parity, spacing, unions and group-duration initialization |
| Compatibility | Source/destination group declarations, categories and ancestry components; handling future templates' categories |
| Calendar inputs | Intensity, template and selection-model schedules; interval columns, coverage, gaps, exact boundaries and clock preservation |
| Households | Party membership, selection probabilities, rate units and reconstruction of internal person/union links |
| Migration history | Stored facts, source provenance, history depth, starting-population unknowns and restart persistence |
| Classifications | Rule interface, output types, update timing, custom-model inputs and sufficient retained information after truncation |
| Rate strata | Declared group versus built-in migration variable, process support and coverage rules |
| Observation | Entry-aware risk/exposure windows, immigrant versus destination-born events, historical support records and population accounting |
| Timing and interface | Arrival tie priority, clock organization, stream identifiers and public constructors |

No implementation milestone or scientific acceptance threshold is established
here. Household implementation can remain deferred while the individual-arrival
contract is developed; its unresolved rate units must not be silently chosen.

## Worked examples for review

These illustrate accepted meanings and proposed validation cases, not executed
results or additional demographic defaults.

**Age and replacement.** A source snapshot ends in 2000 and contains a person
born in 1975. Sampling that person for an arrival in 2050 gives age 25 and
destination birth time 2025. Sampling the same source person again in 2051
gives another age-25 person with birth time 2026 and a different destination
ID. Neither draw removes or ages the source record. If histories accompany
the person, their date translation still needs to be specified consistently.

**Multiple flows.** A stream with intensity 100 and another with intensity
40 have expected individual arrival counts of 100 and 40 in a one-year
interval. Those are expected flows, not quotas or expected final population
increments: arrivals may subsequently die or have children. These flows do
not depend on the destination's current population size.

**Truncation and parity.** At a source endpoint, a living woman has two
recorded children, one alive and one deceased. Retaining living people and
their descendant closure retains both children, so her recorded parity
remains two and her latest recorded birth date remains available. Removing
her deceased ancestors may still require retaining migration-history summaries
and repairing references under the eventual truncation contract.

**Personal status versus family origin.** An immigrant through stream A and
a locally born person with immigrant parents from stream B have a child in
the destination. The child has not personally immigrated and does not inherit
stream A or B as a personal arrival stream. Their family history supports
connections to both streams. A named classification can summarize these
connections for custom inheritance or matching under a user-specified rule;
no universal numeric generation is assigned here.

## Validation plan

Validation must follow the completed specification. For individual arrivals
with a valid, nonempty source, the initial analytical and accounting cases are:

| Case | Expected check |
| --- | --- |
| Constant intensity lambda for T years | Across replications, arrival counts follow Poisson(lambda * T), and waiting times follow Exp(lambda) |
| Zero intensity or zero duration | No arrivals |
| Replacement draws | Repeated source selections get distinct destination IDs; the original template remains unchanged |
| Destination accounting | Final living population equals initial living population plus actual arrivals plus destination births minus destination deaths |
| Imported histories | Historical source births/deaths are excluded from destination event counts and pre-entry destination exposure |
| Truncation | Recorded parity and last-birth dates of retained living women are preserved |

Additional tests should cover multiple streams, age preservation, template
schedule boundaries, future group categories, mixed-stream pedigrees,
classification inputs in custom models, unknown histories, and classification
preservation through truncation and restart. Source and destination should be
reusable across independent seeded runs under the existing RNG contract.
Household checks depend on its deferred sampling and intensity decisions.
Statistical test methods must be justified when implementation begins; no
acceptance thresholds or successful validation results are claimed here.
