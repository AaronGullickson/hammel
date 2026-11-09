# R source organization

Agreed by the user and implemented November 9, 2026 (D143). This guide defines where R functions belong and the rationale for the current structure. Source organization serves maintainers; help topics and website reference groups serve users and are organized independently.

## Placement rules

1. Give each file one responsibility that can be stated in a sentence. Related public functions may share a file; file length is a review signal, not a limit.

2. Keep constructors, convenience setters, and their specific validation together. Keep private helpers near their main consumer unless they implement a clear shared responsibility.

3. Distinguish model specification, execution, observations of population histories, and summaries calculated from input/output rate schedules.

4. Collect observed event-rate calculators in `calculate_rates.R`. Preserve their different event definitions, risk sets, and exposure units. Similar structure does not require a single generic calculator.

5. When a file has a clear lead function, use that function’s exact name as the filename, followed by `.R`. For example, `fertility()` and `set_fertility()` belong in `fertility.R`. Prefer the constructor as the lead when a constructor and convenience setter share a file. For a family of peer functions or shared infrastructure without a natural lead, use a concise name describing the shared responsibility, such as `calculate_rates.R`, `ancestry.R`, or `analysis_utils.R`. Do not add category prefixes or rename functions merely to make filenames sort together. Keep established generated/package metadata filenames.

6. Keep the existing flat `R/` directory and ordinary package namespace. Filenames group source for readers; they do not create isolated modules. Do not introduce manual sourcing or load-order conventions to enforce this map.

7. Keep source organization independent of help topics and pkgdown reference sections. Each demographic process constructor shares a help topic with its convenience setter: `mortality()`/`set_mortality()`, `fertility()`/`set_fertility()`, `union_formation()`/`set_union_formation()`, `union_dissolution()`/`set_union_dissolution()`, and `group_transition()`/`set_group_transition()`. Explain the specification once and distinguish the returned process from the updated model. Keep generic `set_process()` on its own topic, linking to these pairs (D144). The seven observed event-rate calculators share `calculate_rates`, with common arguments documented once and separate sections for their events, risk sets, exposure units, strata, and outputs (D145). Each function keeps its own help alias. Population counts, fertility summaries, and life tables retain separate topics from observed rates. `calculate_tfr()` and `calculate_patfr()` share the `calculate_tfr` topic; `group()` and `set_groups()` share the `group` topic. The two matching-model constructors retain separate topics (D146).

8. Start each hand-maintained R file with a brief comment describing its purpose, above the section headings. Place exported functions first under `# PUBLIC FUNCTIONS ------`, followed by any `# INTERNAL CONSTANTS ------` and `# HELPER FUNCTIONS ------` sections. Omit sections with no contents. Preserve relative order within each group and keep explanatory comments and roxygen blocks with their definitions. Internal-only files need only the applicable internal sections. Keep generated wrappers generated; do not manually reorder or section `RcppExports.R`. Package roxygen metadata uses `# PACKAGE DOCUMENTATION ------`.

9. Use noun-based names for specification constructors: `matching_model()`, `age_gap_matching_model()`, `weight_model()`, and `group_inheritance()`. Use verbs for pipeable actions such as `set_population()`. Retain `create_baseline_population()` as the user-selected exception (D147).

For a new function, first identify its responsibility and closest existing consumer. Extend that file if the responsibility fits. If no existing file is a logical home, explain the proposed responsibility and filename and ask the user before creating a new R script. Do not create one automatically. The user’s authorization of this reorganization covers the files listed below; future additions follow this rule.

## Responsibility boundaries

- **Matching versus union formation:** matching constructors, formula evaluation, and preparation live in `matching_model.R`. Search initiation, eligibility settings, and search-rate validation remain with the process. Group matching retains its own substantial and coherent implementation.

- **Construction versus attachment:** `set_process.R` retains common attachment and dispatch through `attach_process()`. Process-specific constructors and setters live together in their process files. `validate_rates()` remains one shared function, retaining its current name and mortality flag, in `rates_util.R`; splitting its algorithm is separate work.

- **Transition processes versus histories:** process construction, household rules, storage, and resolution stay together. History schema and queries live in `group_history.R`; shared analysis preparation lives in `analysis_utils.R`; union-rate-specific preparation stays with the rate calculators.

- **Observed rates versus demographic summaries:** all seven observed event-rate calculators belong together, including the singularly named `calculate_cohab_convert_rate()`. Keep `calculate_entry_rates()`, `person_rate_rows()`, and `add_unrealized_births()` there as closely related calculation helpers. TFR/PATFR share `calculate_tfr.R`; life tables and population snapshots have distinct files. Shared rate-summary checks live in `summary_rates.R`.

- **Shared analysis helpers:** `analysis_utils.R` contains history preparation, stratification, and output construction shared by rate calculators, `calculate_population_size()`, and `get_detailed_unions()`. Its scope is demographic analysis of recorded histories, including supplied populations. Helpers specific to event-rate calculation, including confidence intervals and exposure-counting wrappers, stay in `calculate_rates.R`.

- **Relationship extraction:** the three small `get_*` relationship functions share `extract_relationships.R`. Kinship remains separate because it has its own traversal contract and shared generation checks. Public names do not change.

- **Cohabitation:** keep its settings/conversion helpers together even though cohabitation is configured through the union-formation process, rather than introducing a new process or public constructor.

- **Simulation:** `simulate()` retains validation, preparation, execution, and result assembly in one function. Its fertility-setting validator lives in `fertility.R` and ancestry encoder in `ancestry.R`. Keep `format_unrealized_births()` with result assembly in `simulate.R`. Extracting preparation/result assembly from the body could be a later refactor after explicit inputs, outputs, validation order, and RNG boundaries have been mapped; it is separate work from file organization.

- **General checks:** retain the three genuinely general checks in `checks.R`. Domain validation stays with its domain; this file must not become a catch-all.

- **Input-rate utilities:** `rates_util.R` combines shared mortality/fertility rate validation, calendar slicing and coverage checks, group-stratum validation, and engine encoding/decoding. Use sections for these related responsibilities. The public proportional-stratification function remains in `stratify_rates.R`; observed-rate calculations remain in `calculate_rates.R`.

- **Database preparation:** `get_database_rates.R` contains both HMD mortality and HFD fertility retrieval/conversion functions and all their import helpers. Use sections for shared request validation and provenance, HMD conversion including open-age pooling, and HFD conversion including parity-specific rates. Neither public importer is the natural lead for both, so use this descriptive filename.

## Current file map

Every top-level R function has one home below. **Bold names are exported**, including `validate_population()`; plain names are internal. Package roxygen metadata stays in `hammel-package.R` (it defines no functions). Shared constants stay with the responsibility they support: `summary_measures` in `summary_rates.R`, `reserved_group_names` in `group.R`, `matching_settings` in `matching_groups.R`, `weight_variable_codes` in `weight_model.R`, and `builtin_strata`/`rate_schema_columns` in `stratify_rates.R`.

### `R/RcppExports.R`

Generated native wrappers; continue to generate rather than hand-edit.

- Functions: `dissolution_counts_cpp()`, `previous_union_info_cpp()`, `previous_union_state_cpp()`, `unpartnered_birth_rows_cpp()`, `union_start_groups_cpp()`, `union_group_history_cpp()`, `extract_kin_cpp()`, `person_lexis_cpp()`, `person_snapshot_cpp()`, `inspect_matching_cpp()`, `inspect_group_matching_cpp()`, `inspect_weighted_calibration_cpp()`, `inspect_exact_draw_cpp()`, `inspect_weighted_incremental_cpp()`, `simulate_cpp()`.

### `R/analysis_utils.R`

Shared history preparation, stratification, and output construction for demographic analysis. Used by rate calculators, population counts, and relationship extraction.

- Functions: `observation_groups()`, `prepare_group_observations()`, `check_by()`, `check_pool_dissolved()`, `previous_labels()`, `person_observation()`, `union_type_labels()`, `build_rows()`, `dim_of()`, `label_dim()`, `person_by_dims()`, `finish_rows()`, `observation_alias()`, `combine_dimensions()`, `finish_strata()`, `rename_group_columns()`.

### `R/ancestry.R`

Ancestry schema, normalization, and initialization.

- Functions: `tracked_ancestry()`, `as_ancestry_list()`, `validate_ancestry()`, `validate_ancestry_table()`, `normalize_population_ancestry()`, `ancestry_matrix()`.

### `R/calculate_life_table.R`

Construct period life tables from mortality rates.

- Functions: **`calculate_life_table()`**.

### `R/calculate_population_size.R`

Observe population counts at specified instants.

- Functions: **`calculate_population_size()`**.

### `R/calculate_rates.R`

Calculate observed event rates from counts and exposure; retain calculator-specific helpers.

- Functions: **`calculate_mortality_rates()`**, `person_rate_rows()`, **`calculate_fertility_rates()`**, `add_unrealized_births()`, `person_breaks()`, `check_demographic_breaks()`, **`calculate_cohab_convert_rate()`**, **`calculate_dissolution_rates()`**, **`calculate_marriage_rates()`**, **`calculate_transition_rates()`**, **`calculate_union_rates()`**, `calculate_entry_rates()`, `prepare_union_group_observations()`, `run_person_lexis()`, `order_rows_keep()`, `band_dim()`, `union_strata()`, `encode_combined()`, `validate_rate_conf_level()`, `add_rate_confidence()`.

### `R/calculate_tfr.R`

Calculate TFR and PATFR from rate schedules.

- Functions: **`calculate_tfr()`**, **`calculate_patfr()`**, `patfr_midpoint()`.

### `R/checks.R`

Small, general checks for table columns and numeric/time values.

- Functions: `check_time()`, `check_columns()`, `check_numeric_column()`.

### `R/cohabitation.R`

Cohabitation and conversion settings within union formation, including preparation.

- Functions: `validate_cohabitation()`, `validate_conversion_rates()`, `prepare_cohabitation()`.

### `R/create_baseline_population.R`

Generate baseline persons and sample declared group categories.

- Functions: **`create_baseline_population()`**, `baseline_group_sampler()`.

### `R/extract_kin.R`

Extract kin and validate generation limits shared with matching.

- Functions: **`extract_kin()`**, `validate_kinship_generations()`.

### `R/extract_relationships.R`

Extract parent, union, and unpartnered-birth relationship tables.

- Functions: **`get_detailed_unions()`**, **`get_parent_info()`**, **`get_unpartnered_births()`**.

### `R/fertility.R`

Construct and attach fertility processes; validate fertility settings and feasibility.

- Functions: **`set_fertility()`**, **`fertility()`**, `unreachable_parity_cells()`, `warn_unreachable_parity_cells()`, `validate_ffab()`, `validate_birth_interval()`.

### `R/get_database_rates.R`

Retrieve HMD mortality and HFD fertility rates, convert them to Hammel tables, and preserve provenance. Include all shared and source-specific import helpers, with sections for shared support, HMD, and HFD.

- Functions: **`get_mortality_hmd()`**, `mortality_from_hmd()`, **`get_fertility_hfd()`**, `hfd_age_intervals()`, `fertility_from_hfd()`, `fertility_by_parity_from_hfd()`, `select_database_year()`, `database_source()`, `check_database_request()`, `check_hmd_shape()`, `hmd_period_intervals()`, `hmd_age_intervals()`, `check_database_years()`, `combine_database_years()`, `check_hmd_open_age()`, `check_hmd_retained_rates()`, `index_hmd_pool_source()`, `pool_hmd_open_age()`.

### `R/group.R`

Declare group dimensions and apply/check declarations.

- Functions: **`group()`**, `validate_group()`, `check_group_names()`, `normalize_group_declaration()`, **`set_groups()`**, `apply_group_declaration()`, `model_groups()`, `declared_group()`, `single_dimension()`, `validate_group_declaration()`.

### `R/group_history.R`

Normalize, validate, and query recorded group histories and dimension views.

- Functions: `empty_transitions()`, `normalize_population_transitions()`, `validate_transition_history()`, `dimension_view()`, `groups_at()`.

### `R/group_inheritance.R`

Construct, validate, and prepare newborn group inheritance rules.

- Functions: **`group_inheritance()`**, `inheritance_methods()`, `normalize_group_inheritance()`, `validate_group_inheritance()`, `extend_identity_groups()`, `prepare_group_inheritance()`.

### `R/group_transition.R`

Construct, attach, validate, and resolve transition processes and household rules.

- Functions: **`group_transition()`**, **`set_group_transition()`**, **`household()`**, `check_transition_group()`, `transition_processes()`, `store_transition_process()`, `resolve_transition_processes()`, `validate_group_transition()`.

### `R/hammel_model.R`

Construct models and attach populations.

- Functions: **`hammel_model()`**, **`set_population()`**, `check_model()`.

### `R/hammel_population.R`

Construct and validate population objects and normalize union tables.

- Functions: `validate_union_types()`, **`hammel_population()`**, **`validate_population()`**, `normalize_unions()`, `empty_unions()`.

### `R/het_fert.R`

Specify persistent fertility heterogeneity and initialize/check individual weights.

- Functions: **`het_fert()`**, `validate_het_fert()`, `fertility_log_variance()`, `validate_fertility_weights()`, `initialize_fertility_weights()`.

### `R/matching_groups.R`

Group pair effects, multidimension matching resolution, and engine representation.

- Functions: `matching_group_labels()`, `expand_matching_pairs()`, `expand_group_exogamy()`, `matching_has_group_formula()`, `matching_predictor_data()`, `matching_pair_offsets()`, `validate_matching_groups()`, `expand_prevent_exogamy()`, `matching_group_allowed()`, `is_dimension_list()`, `split_group_exogamy()`, `exogamy_table_dimensions()`, `validate_exogamy_table()`, `matching_dimensions()`, `internal_matching_model()`, `combined_matching_model()`, `combination_labels()`, `resolve_matching_models()`, `%|na|%()`.

### `R/matching_model.R`

Construct, validate, evaluate, and prepare matching models independently of process settings.

- Functions: **`matching_model()`**, `warn_absorbed_group_formula()`, **`age_gap_matching_model()`**, `validate_delta()`, `validate_matching_model()`, `validate_internal_matching()`, `matching_design_matrix()`, `calculate_matching_utilities()`, `calculate_matching_probabilities()`, `prepare_matching_model()`.

### `R/model_schedules.R`

Validate, select, and prepare calendar schedules of model specifications.

- Functions: `is_model_schedule()`, `validate_model_argument()`, `normalize_model_argument()`, `model_entries()`, `model_at()`, `check_model_coverage()`, `collect_model_bank()`, `prepare_model_schedules()`.

### `R/mortality.R`

Construct and attach mortality processes.

- Functions: **`set_mortality()`**, **`mortality()`**.

### `R/rates_util.R`

Shared input-rate validation, calendar coverage, group strata, and engine encoding/decoding. Keep separate sections within this file for each part of preparing rate schedules for execution.

- Functions: `has_calendar_rates()`, `calendar_rate_segments()`, `calendar_rate_slices()`, `validate_calendar_rates()`, `check_time_coverage()`, `check_calendar_coverage()`, `group_registries()`, `rate_strata_columns()`, `prepare_strata_rates()`, `group_engine_input()`, `strata_spec()`, `decode_report_strata()`, `check_strata_overlap()`, `check_groups()`, `split_rate_groups()`, `validate_rates()`.

### `R/set_process.R`

Attach and revalidate process specifications; construct their common representation.

- Functions: **`set_process()`**, `attach_process()`, `new_rate_process()`.

### `R/simulate.R`

Orchestrate validation, preparation, RNG handling, execution, and result assembly.

- Functions: **`simulate()`**, `format_unrealized_births()`.

### `R/stratify_rates.R`

Public proportional stratification of input rate tables and its checks.

- Functions: **`stratify_rates()`**, `candidate_strata()`, `validate_strata_factors()`, `check_stratified_columns()`.

### `R/summary_rates.R`

Shared strata and interval checks for summaries of rate schedules.

- Functions: `summary_strata()`, `summary_rows()`, `check_summary_intervals()`.

### `R/union_dissolution.R`

Construct and validate dissolution processes.

- Functions: **`union_dissolution()`**, **`set_union_dissolution()`**, `validate_union_dissolution()`.

### `R/union_formation.R`

Construct and validate union formation/search processes.

- Functions: **`union_formation()`**, **`set_union_formation()`**, `validate_union_formation()`, `expand_previous_union_rates()`, `validate_union_search_rates()`.

### `R/weight_model.R`

Construct, validate, and prepare within-cell weight models.

- Functions: **`weight_model()`**, `weight_builtins()`, `weight_variable()`, `weight_group_name_possible()`, `validate_weight_model()`, `weight_design()`, `validate_process_weights()`, `warn_cancelled_weight_terms()`, `weight_strata()`, `prepare_weight_model()`.

## Maintaining the organization

Use the placement rules above when adding functions. Keep constructors and their convenience setters together, and keep domain-specific validation close to the behavior it checks. Ask the user before creating a new R script when no existing file fits.

For a future mechanical reorganization, preserve whole function definitions, constants, explanatory comments, and roxygen blocks. Compare parsed definitions before and after, check package loading, and run `devtools::test()` with live database tests opted out. Regenerate documentation and verify that exports and help topics are preserved. Shared-topic `@order` tags in the mortality/fertility constructors and setters preserve help display order independently of source-file ordering.

Do not combine file moves with function renaming, algorithm changes, defaults, or scientific-contract changes. Any later refactor of `simulate()` needs separate review of its inputs, outputs, validation order, and RNG boundaries, with exact before/after result and RNG comparisons. Preserve user-owned README/vignettes and historical validation reports.

Pkgdown reference grouping remains a separate follow-up choice. The source reorganization preserved help topics; D144 subsequently grouped process constructors with their setters, and D145 combined the observed event-rate help pages as described above. The site’s reference groups are unchanged.

## Verification

The reorganization retains all 226 functions and six shared constants, with unchanged parsed definitions, across 34 R files including generated wrappers and package metadata. All 43 exports are preserved. Original explanatory comments and roxygen text are retained; explicit documentation-order tags keep shared help topics stable. See [testing](testing.md) for verification results.
