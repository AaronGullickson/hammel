# Running tests

The scripts in `tests/` use base R and an installed Hammel package. Run them with `Rscript tests/<name>.R` after installing the package.

## Live database integration

`tests/database-online.R` is opt-in. It downloads USA mortality and fertility rates for 2000, exercises both pipeable database setters, and checks the imported process structure, units, sex schedules, rates, and source metadata. It does not simulate a population or establish demographic accuracy. It never prints credentials or raw database tables.

Install HMDHFDplus and supply `HMD_USERNAME`, `HMD_PASSWORD`, `HFD_USERNAME`, and `HFD_PASSWORD` through your environment or `~/.Renviron`. R reads the user environment file at startup (unless overridden by its startup configuration). Then run:

```sh
HAMMEL_RUN_ONLINE_TESTS=true Rscript tests/database-online.R
```

With no opt-in flag, the test skips before inspecting credentials or loading the optional dependency. With the flag enabled, missing credentials, authentication failures, unavailable data, and conversion errors fail the test. Live data can change; the test checks contracts rather than exact rate values. Offline tests use synthetic fixtures and temporary reader replacements to check exact conversions independently of database availability.

Verified October 12, 2026: authenticated USA/2000 downloads through both setters passed with HMDHFDplus 2.0.8. This verifies the current access and conversion paths, not future database availability or simulation correctness.
