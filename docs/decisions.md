# Decision log

Keep accepted choices separate from proposals. Only record user agreement as acceptance of a demographic or major architectural decision. Each new entry should include status, date, rationale, alternatives, implications, and supporting evidence where available. Update superseded entries with a pointer to their replacement.

## D001: Standalone rewrite

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** build a standalone R package named Hammel, informed by SOCSIM; do not build a SOCSIM wrapper or rsocsim companion package.
- **Rationale:** modernize the user interface and enable custom demographic mechanisms and methodological innovation beyond the legacy file-oriented workflow.
- **Implications:** investigate SOCSIM algorithms and assumptions; legacy input compatibility is optional.
- **Evidence:** project conversation handoff, summarized in `project-brief.md`.

## D002: Scientific ownership

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** the user guides demographic requirements and scientific interpretation; agents support software investigation, implementation, documentation, and validation.
- **Implications:** expose unresolved demographic assumptions and obtain decisions before implementing dependent behavior.

## Pending choices

The design specification lists open choices. In particular, C++, monthly time, continuous time, any example function names, and the four-part candidate architecture are not accepted decisions.
