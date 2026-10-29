# SOCSIM reference audit plan

Updated October 29, 2026. Static source audits cover rates, competing events,
union formation, [dissolution and remarriage](socsim-dissolution-audit.md), and
[cohabitation](socsim-cohabitation-audit.md) for milestone 4. No executable reference validation or comprehensive audit of
all variants has been performed. References below are starting points recorded
in the earlier conversation.

## Establish the reference

An initial static trace of candidate selection, within-month execution, and rescheduling is recorded in [competing events](socsim-competing-events-audit.md).

Record repository URL, commit or release, retrieval date, relevant files, and applicable licenses for each source examined. Distinguish SOCSIM engine code from rsocsim integration. Before copying or adapting code, establish applicable licensing obligations and document provenance. Package licensing remains undecided.

## Reference version and current departures

The existing traces use MPIDR/rsocsim commit
`125450cff968bb3ff762fa99d4c926a8136062a2`. The local checkout was confirmed
at `../rsocsim` relative to the Hammel root during this documentation review;
that location is a convenience, not a portable project dependency.
Source functions and approximate line locations are recorded in the individual
audits. License obligations must be established before copying/adapting code;
Hammel's package license remains undecided and no audited source was copied.

Current Hammel intentionally departs from SOCSIM's monthly timing and input
conventions: continuous chronology, person-year hazards, midpoint fertility
compensation, conditional-logit partner selection, configurable known-kin
screening, and default cancellation of unmatched birth attempts. Source findings do not imply fixed-
seed equivalence or demographic accuracy. Match assumptions/units and use
replicated statistical comparisons when a reference experiment is appropriate.
A [static union-formation trace](socsim-union-formation-audit.md) now covers
queue modes, matching, eligibility, and fertility interactions in the pinned
non-ENHANCED source. It is not executable validation or a complete audit of
all variants.

## Milestone 4 investigation

The user requested an audit of SOCSIM dissolution and remarriage after closing
milestone 3 and creating `milestone4-union-dissolution`. The
[new static trace](socsim-dissolution-audit.md) documents duration-based divorce,
sex/group/status rate indexing, competing spouse candidates, shared-union
closure and partner rescheduling, remarriage status and duration schedules,
implicit rate defaults, cohabitation histories, and former-partner exclusions.
It lists controlled reference experiments that have not yet been run.

After this audit, the user selected duration-based annual hazards with one
female-carried clock per union (D078); see [union dissolution](union-dissolution.md).
SOCSIM still differs in its monthly probability inputs and ability to give both
spouses a risk. Its implicit rate defaults, reconciliation rules and queue
quirks are not automatically adopted. The selected cohabitation behavior is
specified separately under D085–D086.

## Cohabitation stretch-goal investigation

On October 28, after the user described milestone 4's required elements as
wrapped up, a focused [cohabitation audit](socsim-cohabitation-audit.md) traced
formation/type probabilities, same-partner legalization, age/duration conversion
schedules, breakup/death and earlier legal-status recovery, fertility, implicit
defaults, record types and counters in the same pinned source. Static caveats
include a missing-looking one-queue male re-entry after cohabitant death and
mixed/omitted legal-marriage event counters. Those need controlled reference
runs before their population-level effects can be claimed. No reference
simulation was executed, no code copied, and no Hammel cohabitation design
was selected by the audit itself. Subsequent D085–D086 specify and implement
the rate-mode design, and the user accepted milestone 4 on October 29, 2026
(D093). Static reference findings and suggested experiments remain distinct
from the accepted Hammel validation report.

## Parity-specific birth-spacing investigation

The [parity/spacing trace](socsim-parity-spacing-audit.md), inspected October 29,
2026 at the same pinned commit, verifies SOCSIM's previous-parity adjustment,
unadjusted parity zero, actual history-derived parity and highest-parity copy.
It is a static trace, not executable evidence or an accepted Hammel correction.

## Reference links

Starting references:

- SOCSIM/rsocsim: https://github.com/MPIDR/rsocsim
- Source tree: https://github.com/MPIDR/rsocsim/tree/main/src
- R interface documentation: https://mpidr.github.io/rsocsim/reference/socsim.html
- Workshop materials: https://github.com/alburezg/rsocsim_workshop_paa
- Berkeley SOCSIM site: https://lab.demog.berkeley.edu/socsim/

## Audit areas

| Area | Questions to answer |
| --- | --- |
| Scheduling | How are events sampled, ordered, competed, canceled, and rescheduled? |
| Rates | What do inputs mean, in what units, and how are age and time intervals handled? |
| Mortality | Who is at risk, how is death scheduled, and what dependent state changes follow? |
| Fertility | What determines eligibility, heterogeneity, birth timing, and parent assignment? |
| Partnerships | How are candidates selected, matched, and constrained; how are unions dissolved? |
| Kinship | How are identities, parenthood, relationship histories, and consistency maintained? |
| State and RNG | What is per-run state, what is global, how are seeds handled, and what survives between runs? |
| Initialization and output | What histories or burn-in are required, and when are observations recorded? |

## Finding template

For each finding, record:

- Pinned source version and file/function/line locations.
- Behavior in demographic terms, with a worked example where useful.
- Evidence: source trace, documentation, or an executable reference experiment.
- Uncertainties and gaps, including disagreement between code and documentation.
- Consequences for a candidate Hammel design; mark these as proposals.
- Validation cases suggested by the finding.

Do not treat legacy behavior as automatically desirable or scientifically correct. Separate confirmed behavior from interpretation. The earlier report that in-process runs can affect later simulations is an investigation lead, not a diagnosed cause.
