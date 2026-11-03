# Specify persistent fertility weights and maternal inheritance

Weights redistribute births among eligible women within a fertility
cell; they do not multiply the cell's target rate. Fresh weights are
lognormal with mean one and coefficient of variation `cv`. Daughters
inherit on the normal scale: Z_d = alpha \* Z_m + sqrt(1 - alpha^2) \*
epsilon. Here alpha describes log-weight association, not raw-weight
correlation. This preserves the lognormal distribution when reproducing
mothers have the specified latent normal distribution; reproductive
selection can change it across generations.

## Usage

``` r
het_fert(cv, alpha)
```

## Arguments

- cv:

  Finite nonnegative scalar coefficient of variation. Zero means fresh
  weights are one, with no random draw.

- alpha:

  Finite scalar in `[0,1]`. Zero gives independent daughters; one copies
  maternal weights exactly. At zero CV the log-space formula gives W_d =
  W_m^alpha for supplied nonunit weights, without random innovation.

## Value

A `hammel_het_fert` specification for
[`set_fertility()`](https://aarongullickson.github.io/hammel/reference/set_process.md)
or
[`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md).
Construction consumes no random draws.

## Details

Simulation draws weights for all recorded women when persons has no
`hetfert` column, or retains supplied finite strictly positive female
values. An incomplete supplied column is an error. Male weights are
unused; newly created male entries are NA. Draws follow the simulation
seed. The output persons table always contains `hetfert` when this model
is enabled, including zero-year runs. Enabled restarts retain weights.
With `hetfert = NULL`, output omits the column even if supplied, so
later re-enabling draws fresh weights for all women.

## Examples

``` r
het_fert(cv = 0.5, alpha = 0.4)
#> $cv
#> [1] 0.5
#> 
#> $alpha
#> [1] 0.4
#> 
#> attr(,"class")
#> [1] "hammel_het_fert"
```
