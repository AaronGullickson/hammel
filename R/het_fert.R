# Specify persistent fertility heterogeneity and initialize individual weights.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify persistent fertility weights and maternal inheritance
#'
#' @description
#'
#' Use fertility weights that allow for greater inheritable variance in
#' fertility behavior across women.
#'
#' @details
#'
#' With het_fert in use, a `hetfert` column is added to the persons table with a
#' value for every woman. This value \eqn{Z} is drawn from a standard normal
#' distribution. That value is converted into a weight \eqn{W}:
#'
#' \deqn{
#'   W = e^{-\sigma^2 / 2 + \sigma Z},\qquad\sigma^2 = log(1 + c^2)
#' }
#'
#' where \eqn{c} is the coefficient of variation. \eqn{W} will have a lognormal
#' distribution with \eqn{E(W)=1} and \eqn{SD(W)=c}. The weight \eqn{W} is
#' applied to the selection of eligible woman for a birth event within a
#' fertility rate cell.
#'
#' Daughters may partially inherit their score \eqn{Z_d}. The daughter's latent
#' normal score \eqn{Z_d} is determined by her own draw from the standard normal
#' distribution \eqn{Z_a} and the mother's score as follows:
#'
#' \deqn{
#'   Z_d = \alpha Z_m + \sqrt{1-\alpha^2}Z_a
#' }
#'
#' where \eqn{\alpha} determines the degree of inheritance. A value of zero
#' means no heritability while a value of 1 means perfect heritability.
#'
#' The correlation coefficient between mother's and daughter's weights is given
#' by:
#'
#' \deqn{
#'  r(W_m, W_d)=\frac{(1+c^2)^\alpha-1}{c^2}
#' }
#'
#' @param cv Finite nonnegative scalar coefficient of variation that defines
#'   amount of variation in fertility. Zero indicates no variance.
#' @param alpha Finite scalar in `[0,1]` that defines degree of inheritability.
#'   Zero gives independent daughters; one copies maternal weights exactly.
#' @return A `hammel_het_fert` specification for [set_fertility()] or
#'   [fertility()].
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
  if (is.null(model)) {
    return(invisible(NULL))
  }
  if (!inherits(model, "hammel_het_fert") || !is.list(model)) {
    stop(
      "`hetfert` must be NULL or a het_fert(cv, alpha) object.",
      call. = FALSE
    )
  }
  scalar <- function(x) {
    is.numeric(x) && is.null(dim(x)) && length(x) == 1L &&
      is.finite(x)
  }
  if (!scalar(model$cv) || model$cv < 0) {
    stop("`cv` must be a finite nonnegative numeric scalar.", call. = FALSE)
  }
  if (!scalar(model$alpha) || model$alpha < 0 || model$alpha > 1) {
    stop("`alpha` must be a finite numeric scalar in [0, 1].", call. = FALSE)
  }
  invisible(model)
}

fertility_log_variance <- function(cv) {
  if (cv <= 1) log1p(cv^2) else 2 * log(cv) + log1p((1 / cv)^2)
}

validate_fertility_weights <- function(persons) {
  if (!"hetfert" %in% names(persons)) {
    return(invisible(NULL))
  }
  if (!"sex" %in% names(persons) || anyNA(persons$sex)) {
    stop(
      paste0(
        "Fertility weight validation requires nonmissing ",
        "`persons$sex` values."
      ),
      call. = FALSE
    )
  }
  w <- persons$hetfert
  # Base R commonly represents an all-NA column as logical. It is a valid
  # unused placeholder in male-only/empty histories, but not for women.
  if (!is.null(dim(w)) ||
        !(is.numeric(w) || (is.logical(w) && all(is.na(w))))) {
    stop(
      "`persons$hetfert` must be a numeric column; male entries may be NA.",
      call. = FALSE
    )
  }
  female <- which(persons$sex == "f")
  if (any(!is.finite(w[female]) | w[female] <= 0)) {
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
  }
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
    persons$hetfert[female] <- if (model$cv == 0) {
      rep(1, length(female))
    } else {
      stats::rlnorm(
        length(female),
        meanlog = -variance / 2,
        sdlog = sqrt(variance)
      )
    }
  }
  validate_fertility_weights(persons)
  persons$hetfert <- as.numeric(persons$hetfert)
  persons
}
