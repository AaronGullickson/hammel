# Hammel project brief

Status: project direction recorded from the October 11, 2026 planning conversation. The initial package skeleton exists; simulation behavior has not been implemented or audited.

## Accepted direction

Hammel, named in honor of one of SOCSIM's creators, will be a standalone demographic microsimulator in an R package. SOCSIM provides a reference implementation and scientific experience to investigate and adapt. Hammel is a substantial rewrite, not a wrapper or companion package for rsocsim.

The intended user experience includes:

- R functions to construct populations, demographic inputs, and processes.
- Pipeable model construction and simulation workflows.
- R-native results that users can inspect and analyze.
- Custom demographic mechanisms and opportunities for methodological innovation.
- A maintainable implementation whose demographic assumptions can be examined and validated.

The user owns demographic requirements and scientific interpretation. Agents support source investigation, implementation, documentation, and validation. Research use requires evidence beyond successful execution or plausible populations.

## Boundaries

Legacy SOCSIM file compatibility is optional and must not govern the new design. A universal simulation framework is not a project requirement. The implementation language, dependencies, precise API, and demographic algorithms remain undecided.

## First milestone

Accepted October 12, 2026: simulate a population with demographic rates that remain constant over calendar time, then inspect results and test demographic accuracy. Fertility and mortality are the focus. Include a basic marriage or partnership mechanism to assign parents; a realistic partnership market is outside this milestone.

Validate fertility and mortality against the supplied inputs, including total fertility rate (TFR) and life-table measures. Constant over calendar time does not imply constant across ages; the age pattern and rate representation still need specification. Random pairing after a minimum age is a candidate mechanism, not yet an accepted algorithm.

## Proposed development sequence

These milestones are planning suggestions, not implementation authorization:

1. Establish reference source versions and licensing; audit relevant SOCSIM behavior.
2. Resolve enough of the design and scientific contract to specify a minimal simulator.
3. Build the first milestone incrementally: supplied-population initialization and mortality, then fertility, births, parent–child relationships, and basic partnerships.
4. Validate output fertility and life-table measures against the inputs, alongside accounting and relationship checks.
5. Extend the validated foundation with richer partnership behavior and transitions.

Keep specifications, worked examples, decisions, audit evidence, and validation results in this repository so subsequent agents can resume from durable evidence.

## Provenance

This brief and the initial companion documents were distilled from `/home/aarong/Work/Hammel-project-notes.md`, saved October 11, 2026. That file summarizes a planning conversation; it does not contain a detailed source audit. These repository documents are the working guidance going forward. The decision log records later changes.

## Monthly time decision

The first milestone uses a fixed monthly simulation time unit, following SOCSIM's convention. Configurable units are a longer-term goal. Age and event-date conventions, exposure accounting, within-month ordering, and conversion of input rates to monthly event probabilities or schedules still require specification. A monthly time unit does not imply simply dividing annual probabilities by twelve.

## User-facing rate units

Users can supply age-specific fertility and mortality rates in events per person-year. The monthly engine converts these to events per person-month internally by dividing by twelve. This removes the need for users to prepare monthly rates while preserving explicit rate interpretation and denominators.
