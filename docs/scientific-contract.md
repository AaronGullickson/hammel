# Scientific contract

Status: requirements for specifying and validating behavior. No mortality, fertility, partnership, or scheduling algorithm has yet been accepted.

## General requirements

Every consequential mechanism needs a plain-language definition and worked examples before implementation. Explicitly distinguish hazards, probabilities, and central rates. State units, interval definitions, conversions, and interpolation assumptions. Do not apply a conversion merely because it is conventional: identify the assumptions that make it valid.

Specifications must define eligibility, event timing, competing-event semantics, state transitions, and consequences for other pending events. Record intentional differences from SOCSIM and their scientific implications.

## Mechanism specification template

Use this structure for each proposed mechanism:

1. **Status and scope:** proposed or accepted; behavior covered and excluded.
2. **Inputs:** schema, units, valid ranges, missing-value policy, and provenance.
3. **Eligibility:** population at risk and dependence on age, history, relationships, or attributes.
4. **Timing:** rate or probability interpretation, scheduling, interval boundaries, ties, and competition.
5. **Transition:** person and relationship changes; dependent events invalidated or rescheduled.
6. **Observations:** emitted events, recorded attributes, and snapshot consequences.
7. **Worked examples:** ordinary cases and boundary cases with expected outcomes.
8. **Validation:** analytical expectations, accounting and relationship invariants, and reference comparisons.
9. **Open questions:** decisions needed from the user before dependent implementation.

## Foundation requiring specification

The proposed first milestone is supplied-population initialization and mortality. Before implementation, resolve input requirements, time semantics, mortality input interpretation, competing events, observation timing, and RNG behavior. Fertility then requires explicit eligibility, heterogeneity, newborn initialization, and parenthood rules. Partnerships require matching and market behavior specifications.

## Validation principles

- Derive checks from the scientific specification, rather than reproducing implementation logic.
- Use analytical cases with known answers, including deterministic boundary cases where applicable.
- Check population accounting: for a closed population, living population changes equal births minus deaths. Extend the identity explicitly when migration or other entry and exit mechanisms are introduced.
- Check valid event dates, person identities, and relationship consistency under the chosen model.
- Verify reproducibility under the declared seed contract and independence between simulations.
- Compare demographic outcomes against a pinned SOCSIM reference over multiple replications when relevant. Define comparison measures and acceptance criteria before interpreting results.
- Require exact fixed-seed comparisons only when algorithms and random-number consumption are intentionally preserved. For changed algorithms, use statistical comparisons with justified uncertainty and tolerances.

Passing software checks does not establish research readiness. Before research use, seek targeted review of scheduling, RNG handling, and the R/native boundary by an experienced simulation developer, alongside demographic review.
