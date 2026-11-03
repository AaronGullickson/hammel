# Specify individual group transitions

Construct or attach competing origin-to-destination hazards. Only the
individual changes group; ancestry, parentage and union histories are
retained. Duration is years since the latest recorded transition, or age
without one. Missing directed pairs, strata and gaps have zero risk.
Restarts retain transition histories and redraw pending events.

## Usage

``` r
group_transition(rates)

set_group_transition(model, rates)
```

## Arguments

- rates:

  A nonempty data frame requiring `group_origin`, `group_dest`, `rate`,
  and either `age_start`/`age_end` or `dur_start`/`dur_end` throughout.
  Half-open bands use years; rates are finite nonnegative hazards per
  person-year. Ends may be `Inf`. Bands cannot overlap within a stratum.
  Self-transitions are invalid. Optional `sex` is m/f. Optional
  `union_type` is unpartnered/partnered or unpartnered/cohab/marriage.
  Omitted columns apply universally. Destinations may initially be
  empty; other configured group-specific processes must cover all
  possible groups. Optional `time_start`/`time_end` specify half-open
  absolute calendar-year bands, including fractional years. Every active
  slice must be valid and the table must cover the requested run. Rate
  changes do not reset group duration.

- model:

  A `hammel_model`.

## Value

`group_transition()` returns a process with a rates tibble;
`set_group_transition()` returns the updated model. Grouped simulation
output always includes a `transitions` tibble with `person_id`, `time`,
`group_origin`, and `group_dest`, including when empty.
