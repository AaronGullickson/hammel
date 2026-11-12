# Create a weight model for selection within rate cells

Create a `weight_model` that will be used to generate weights for
determining which individuals/units experience an event.

## Usage

``` r
weight_model(formula = NULL, coefficients = numeric(), fun = NULL)
```

## Arguments

- formula:

  A one-sided formula of weight variables.

- coefficients:

  A finite named numeric vector of log hazard ratios.

- fun:

  A custom function returning log weights, instead of `formula`.

## Value

A `hammel_weight_model`. This specifies a model; it does not estimate
one.

## Details

A weight model decides which member of a rate cell experiences an event,
without changing the cell's input intensity and thus the input rates.
The cell produces events at its input rate times its members; each event
goes to member `i` with probability `w_i / sum(w)`, where
`w_i = exp(eta_i)`.

With a formula, `eta_i` is the linear predictor `x_i %*% beta`, so
`exp(beta)` is the hazard ratio between two members of the same cell
that differ by one unit in a term. Users can therefore create
proportional hazards models, without affecting overall input rates.
Population level effects can be smaller where cells are small or
homogeneous.

Attach a weight model with the `weights` argument of
[`set_mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
[`set_fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
[`set_union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md),
[`set_group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)
and
[`set_union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md),
or `weights_conversion` of
[`set_union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md).

## Variables

Group variables are named by the declared group dimensions (see
[`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md));
for a dimension `race` with categories A and B.

Mortality, fertility and group transitions are person processes and may
use `sex` (`"f"`, `"m"`), `race`, `union_type` (`"unpartnered"`,
`"cohab"`, `"marriage"`), `partnered` (0/1), `parity` (births; 0 for
men), `previous_union` (`"none"`, `"separated"`, `"widowed"`),
`partner_race` (NA when unpartnered) and, when the dimension tracks
ancestry, shares `ancestry_race_A` and `ancestry_race_unknown`. Group
transitions may also use `destination`, the destination of the risk in
the transition's dimension.

Dissolution and conversion are union processes and may use `race_w`,
`race_h`, `parity_w`, `previous_union_w`, `previous_union_h`,
`ancestry_race_A_w`, `ancestry_race_A_h`, `same_race` (0/1),
`age_difference` (husband's age minus wife's age) and `union_type`
(`"cohab"`, `"marriage"`).

The first listed level of each categorical variable is its reference;
for a group dimension it is its first category (declared, or first in
the population).

## Formulas

Terms are main effects of these variables and interactions (`:` or `*`);
transformations such as [`I()`](https://rdrr.io/r/base/AsIs.html) are
not supported (use `fun`). Categorical variables use treatment coding,
and `coefficients` must be named exactly as the corresponding
model-matrix columns, for example `raceB` or `union_typecohab`.

## Custom functions

Instead of using a formula, users can create a custom function with
`fun`. `fun(members, bin)` receives a tibble with one row per eligible
member of the chosen cell, holding the variables above for every
declared dimension (with labels for categorical variables) and
identifiers (`person_id`, or `union_id`, `wife_id` and `husband_id`),
and a list `bin` describing the cell: its `process`, strata (each
dimension's category, or `<dim>_w`/`<dim>_h` for union processes, NA
when the rates do not stratify by it; `origin` and `destination` for
transitions), band (`age_start`/`age_end` or `dur_start`/`dur_end`) and
the active rate period `time_start`/`time_end` (NA for static rates).
Exact ages and times are not passed. In fertility with weighted matching
calibration, the function is evaluated at every refresh of every woman,
which can be slow.

## Time-varying weights

A process's `weights` may instead be a tibble with numeric `time_start`,
`time_end` and a list column `model` of weight models, as for matching
and inheritance model schedules.

## Examples

``` r
# allow mortality to be 10% lower for individuals with partners
partnered_mortality <- weight_model(~partnered, c(partnered = log(0.9)))

# allow dissolution rates to change by 5% per year increase in partner
# age differences - not in absolute terms which would require the custom
# function
gap_dissolution <- weight_model(
  ~age_difference,
  c(age_difference = log(1.05))
)
# like previous example but use custom function to get absolute age
# difference
gap_dissolution <- weight_model(
  fun = function(members, bin) {
    log(1.05) * abs(members$age_difference)
  }
)
```
