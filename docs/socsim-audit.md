# SOCSIM reference audit plan

Status: planned. No detailed source audit has been performed. References below are starting points recorded in the earlier conversation, not verified findings in this repository.

## Establish the reference

Record repository URL, commit or release, retrieval date, relevant files, and applicable licenses for each source examined. Distinguish SOCSIM engine code from rsocsim integration. Before copying or adapting code, establish applicable licensing obligations and document provenance. Package licensing remains undecided.

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
