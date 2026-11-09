# Specify persistent fertility heterogeneity and initialize individual weights.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify persistent fertility weights and maternal inheritance
#'
#' Weights redistribute births among eligible women within a fertility cell;
#' they do not multiply the cell's target rate. Fresh weights are lognormal with
#' mean one and coefficient of variation `cv`. Daughters inherit on the normal
#' scale: Z_d = alpha * Z_m + sqrt(1 - alpha^2) * epsilon. Here alpha describes
#' log-weight association, not raw-weight correlation. This preserves the
#' lognormal distribution when reproducing mothers have the specified latent
#' normal distribution; reproductive selection can change it across generations.
#'
#' Simulation draws weights for all recorded women when persons has no `hetfert`
#' column, or retains supplied finite strictly positive female values. An
#' incomplete supplied column is an error. Male weights are unused; newly
#' created male entries are NA. Draws follow the simulation seed. The output
#' persons table always contains `hetfert` when this model is enabled, including
#' zero-year runs. Enabled restarts retain weights. With `hetfert = NULL`,
#' output omits the column even if supplied, so later re-enabling draws fresh
#' weights for all women.
#'
#' @param cv Finite nonnegative scalar coefficient of variation. Zero means
#'   fresh weights are one, with no random draw.
#' @param alpha Finite scalar in `[0,1]`. Zero gives independent daughters; one
#'   copies maternal weights exactly. At zero CV the log-space formula gives W_d
#'   = W_m^alpha for supplied nonunit weights, without random innovation.
#' @return A `hammel_het_fert` specification for [set_fertility()] or
#'   [fertility()]. Construction consumes no random draws.
#' @export
#' @examples
#' het_fert(cv = 0.5, alpha = 0.4)
het_fert <- function(cv, alpha) {
  model <- structure(list(cv = cv, alpha = alpha), class = "hammel_het_fert")
  validate_het_fert(model)
  model
}

# HELPER FUNCTIONS ---------------------------------------------------------

validate_het_fert <- function(model) {
  if (is.null(model)) return(invisible(NULL))
  if (!inherits(model, "hammel_het_fert") || !is.list(model))
    stop("`hetfert` must be NULL or a het_fert(cv, alpha) object.",
      call. =
        FALSE
    )
  scalar <- function(x) is.numeric(x) && is.null(dim(x)) && length(x) == 1L &&
    is.finite(x)
  if (!scalar(model$cv) || model$cv < 0)
    stop("`cv` must be a finite nonnegative numeric scalar.", call. = FALSE)
  if (!scalar(model$alpha) || model$alpha < 0 || model$alpha > 1)
    stop("`alpha` must be a finite numeric scalar in [0, 1].", call. = FALSE)
  invisible(model)
}

fertility_log_variance <- function(cv) {
  if (cv <= 1) log1p(cv^2) else 2 * log(cv) + log1p((1 / cv)^2)
}

validate_fertility_weights <- function(persons) {
  if (!"hetfert" %in% names(persons)) return(invisible(NULL))
  if (!"sex" %in% names(persons) || anyNA(persons$sex))
    stop(
      paste0(
        "Fertility weight validation requires nonmissing ",
        "`persons$sex` values."
      ),
      call. = FALSE
    )
  w <- persons$hetfert
  # Base R commonly represents an all-NA column as logical. It is a valid
  # unused placeholder in male-only/empty histories, but not for women.
  if (!is.null(dim(w)) || !(is.numeric(w) || (is.logical(w) && all(is.na(w)))))
    stop("`persons$hetfert` must be a numeric column; male entries may be NA.",
      call. = FALSE
    )
  female <- which(persons$sex == "f")
  if (any(!is.finite(w[female]) | w[female] <= 0))
    stop(
      paste(
        paste0(
          "Female `persons$hetfert` values must be finite and ",
          "strictly positive."
        ),
        "Remove the column to draw fresh weights for all women."
      ),
      call. = FALSE
    )
  invisible(NULL)
}

initialize_fertility_weights <- function(persons, model) {
  if (is.null(model)) {
    persons$hetfert <- NULL
    return(persons)
  }
  if (!"hetfert" %in% names(persons)) {
    female <- which(persons$sex == "f")
    persons$hetfert <- rep(NA_real_, nrow(persons))
    variance <- fertility_log_variance(model$cv)
    persons$hetfert[female] <- if (model$cv == 0) rep(1, length(female)) else
      stats::rlnorm(length(female), meanlog = -variance / 2, sdlog = sqrt(
        variance
      ))
  }
  validate_fertility_weights(persons)
  persons$hetfert <- as.numeric(persons$hetfert)
  persons
}
