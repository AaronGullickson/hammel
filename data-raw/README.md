# Data inputs and preparation

Data, original source files, provenance and offline preparation scripts are
grouped by country and period:

| Folder | Contents |
| --- | --- |
| [england_wales_1992_onward](england_wales_1992_onward/README.md) | Partnership, fertility and divorce evidence for the second-transition demonstration |
| [usa_2000_2001](usa_2000_2001/usa_2000_rates.README.md) | Fertility in 2000, marriage transitions/rates in 2001, original sources and preparation script; ballpark [regional migration rates](usa_2000_2001/usa_2000_regional_migration.README.md) by age, around 2000; [1999-2001 mortality](usa_2000_2001/usa_1999_2001_mortality.README.md) abridged from NCHS decennial life tables. The package's `usa_2000_*` datasets are built from these files by `data-raw/usa_2000.R` |
| [niger_2012](niger_2012) | WPP 2024 age-specific fertility estimates for Niger in 2012, with source metadata |

Preparation scripts run from the package root and now read/write within their
respective folders. Individual filenames and numerical data are unchanged.

Historical validation reports use the input paths from their original revisions.
