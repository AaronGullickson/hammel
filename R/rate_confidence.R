# Shared by public event-rate summaries, after their final aggregation.
validate_rate_conf_level <- function(conf.level) {
  if (!is.null(conf.level) && (!is.numeric(conf.level) ||
      length(conf.level) != 1L || !is.null(dim(conf.level)) ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1))
    stop("`conf.level` must be NULL or a finite numeric scalar strictly between 0 and 1.",
      call. = FALSE)
  invisible(NULL)
}

add_rate_confidence <- function(table, events, exposure, conf.level) {
  validate_rate_conf_level(conf.level)
  if (is.null(conf.level)) return(table)
  table$conf_low <- rep(NA_real_, nrow(table))
  table$conf_high <- rep(NA_real_, nrow(table))
  for (i in which(table[[exposure]] > 0)) {
    limits <- stats::poisson.test(x = table[[events]][i], T = table[[exposure]][i],
      conf.level = conf.level)$conf.int
    table$conf_low[i] <- limits[1L]
    table$conf_high[i] <- limits[2L]
  }
  table
}
