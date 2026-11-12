# Download the HMD and HFD rates used by the website articles and save them
# to vignettes/articles/data, so the site can build without network access
# or database credentials. Run from the package root with hammel installed and
# HMD_USERNAME, HMD_PASSWORD, HFD_USERNAME and HFD_PASSWORD set (see
# ?get_mortality_hmd). Rerun to refresh the saved rates.
#
# Sources and license: Human Mortality Database (https://www.mortality.org/)
# and Human Fertility Database (https://www.humanfertility.org/), both
# licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
# The saved tables are derived from those data: get_mortality_hmd() estimates
# the open age group from pooled deaths and exposures, and both functions
# reshape the series into hammel rate schedules. Each table keeps its download
# details, including the access date, in attr(x, "source"). The articles make
# further changes after loading; those are shown in the articles.

library(hammel)

out_dir <- "vignettes/articles/data"
dir.create(out_dir, showWarnings = FALSE)

save_rates <- function(rates, file) {
  stopifnot(nrow(rates) > 0L, !is.null(attr(rates, "source")))
  saveRDS(rates, file.path(out_dir, file))
}

# demographic_transition.qmd
save_rates(
  get_fertility_hfd("SWE", 1891:2024),
  "sweden_fertility.rds"
)
save_rates(
  get_mortality_hmd("SWE", 1751:2024, open_age = 85),
  "sweden_mortality.rds"
)

# second_demog_transition.qmd
save_rates(
  get_mortality_hmd("GBRTENW", 1992, shape = "5x1", open_age = 90),
  "england_wales_mortality_1992.rds"
)
