# Niger fertility, 2012

`niger_fertility` contains single-year age-specific fertility estimates for
Niger in 2012 from the United Nations World Population Prospects 2024 revision.
The retained CSV covers ages 15–49 and expresses rates in births per woman-year.
The source rates were divided by 1,000. The sum over those ages is 7.468122
births per woman.

[sources.json](sources.json) records the original download, selection,
conversion, and checksums for the download and retained extract.

Run the builder from the package root with `tibble` installed.

```sh
Rscript data-raw/niger_2012/niger_fertility.R
```
