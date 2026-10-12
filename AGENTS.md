# Working on Hammel

Read [the project brief](docs/project-brief.md), [the design specification](docs/design.md), and [the scientific contract](docs/scientific-contract.md) before consequential implementation. Consult [the decision log](docs/decisions.md) for the status and rationale of choices and [the audit plan](docs/socsim-audit.md) when investigating SOCSIM.

## Purpose and scope

Hammel is a new, standalone demographic microsimulator with a modern R package interface, informed by SOCSIM's approach and features, but offering more flexibility to users to create customizable model parameters. It is a substantial rewrite, not a SOCSIM wrapper, rsocsim companion package, or upstream contribution. Legacy file compatibility must not determine its architecture.

The user guides demographic requirements and scientific interpretation. Agents investigate, implement, document, and validate software. Do not silently resolve questions that change the demographic model.

## Before implementing a feature

1. Identify the requested behavior and the relevant accepted decisions.
2. For consequential simulation behavior, write or update its plain-language scientific specification and worked examples first. Identify units, eligibility, event timing, edge cases, and expected outcomes.
3. Identify unresolved modeling choices. Explain their consequences to the user and obtain a decision before implementing behavior that depends on them. Continue independent investigation or infrastructure work where possible.
4. Define validation from the specification, including analytical cases or invariants where applicable.
5. Implement a narrow, reviewable change and keep the documents consistent with the resulting behavior.

Suggestions, illustrative APIs, and proposed milestones are not accepted decisions or authorization to implement them. Routine software choices within an authorized task can be made autonomously; record choices with substantial architectural consequences.

## Engineering expectations

- Keep model construction separate from simulation execution. Aim for composable R workflows and R-native, inspectable results at the user-facing boundary. This does not require an R-only implementation: performance-critical computation and internal state may use compiled code, such as C++, with results exposed as ordinary R objects for inspection and analysis.
- Each simulation must own independent state. Avoid hidden process-wide population state and undocumented effects on later runs.
- Make time, rate units, random-number behavior, and observation semantics explicit.
- Design built-in and custom demographic mechanisms around a documented contract; arbitrary callbacks must not bypass state consistency or event invalidation rules.
- Preserve scientifically relevant provenance. Before copying or adapting source, establish its version, license, and applicable obligations.
- Do not select a compiled implementation language or add major dependencies merely because the notes mention a candidate.
- Do not invent demographic defaults, fabricate audit findings, or claim scientific validity from plausible-looking output.
- Do not start sub-agents solely because the project notes describe possible agent roles. Use delegation when the current task authorizes it.

## Validation and reporting

Use tests appropriate to the behavior: analytical cases, accounting identities, valid dates and relationships, reproducibility, and independence between runs. Compare with SOCSIM over multiple replications when a reference comparison is appropriate. Exact fixed-seed agreement is meaningful only when algorithms and random-number consumption are preserved; otherwise use justified statistical comparisons.

Report what changed, what was checked, and remaining limitations. Distinguish implemented behavior, verified findings, proposals, and unresolved assumptions. Document intentional departures from SOCSIM. Never mark a scientific decision accepted without the user's agreement.
