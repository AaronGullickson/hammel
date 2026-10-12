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

## Proposed development sequence

These milestones are planning suggestions, not implementation authorization:

1. Establish reference source versions and licensing; audit relevant SOCSIM behavior.
2. Resolve enough of the design and scientific contract to specify a minimal simulator.
3. Implement supplied-population initialization and mortality through the intended R interface.
4. Add births and parent–child relationships, exercising the process extension contract.
5. Validate that foundation before adding partnerships and richer transitions.

Keep specifications, worked examples, decisions, audit evidence, and validation results in this repository so subsequent agents can resume from durable evidence.

## Provenance

This brief and the initial companion documents were distilled from `/home/aarong/Work/Hammel-project-notes.md`, saved October 11, 2026. That file summarizes a planning conversation; it does not contain a detailed source audit. These repository documents are the working guidance going forward. The decision log records later changes.
