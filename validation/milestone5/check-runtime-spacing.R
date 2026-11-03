# Offline method-investigation diagnostic, not Hammel engine code or an
# acceptance test. Independent discrete-time toy comparing spacing corrections
# when the population mix differs from the one assumed in advance.
# Run from the repository root:
#   Rscript validation/milestone5/check-runtime-spacing.R
# Writes validation/milestone5/runtime-spacing-summary.csv.
#
# Initialized cohort of 40,000 women at parity zero, spacing B = .75, no
# mortality, one replication (seed 7). Conditional targets per woman-year at
# current parity are .2/.3/.4 (2+) for t < 5 and 1.6 times those for t >= 5.
# Methods:
#   none        : control with no spacing restriction at all
#   composition : h_k = m_k / e_k(t), e_k from the deterministic target cohort
#                 built from the ORIGINAL schedule (see check-composition-spacing.R)
#   runtime     : cell = input parity band; each eligible woman's hazard is
#                 m_c * N_c / N_eligible_c, using current simulated counts.
#                 Equivalent to cell intensity m_c * N_c assigned to a random
#                 eligible woman. Intensity in cells with no eligible women is
#                 recorded as expected lost births.
# See design/parity-spacing-investigation.md for interpretation and limitations.
set.seed(7)
B <- .75; dt <- .0025; Tend <- 8; K <- 15; N <- 40000
base <- function(k) ifelse(k == 0, .2, ifelse(k == 1, .3, .4))
mult <- function(t) ifelse(t < 5, 1, 1.6)
steps <- round(Tend / dt); tt <- (0:steps) * dt
group <- function(p) pmin(p, 2L) + 1L

mk <- base(0:K); l <- matrix(0, steps + 1, K + 1); l[1, 1] <- 1
for (i in 1:steps) {
  out <- mk * l[i, ]; out[K + 1] <- 0
  l[i + 1, ] <- l[i, ] - dt * out + dt * c(0, out[-(K + 1)])
}
inflow <- cbind(0, sweep(l[, -(K + 1)], 2, mk[-(K + 1)], `*`))
cumin <- apply(inflow, 2, cumsum) * dt
lagi <- pmax(0, (0:steps) - round(B / dt)) + 1
elig <- 1 - (cumin - cumin[lagi, ]) / pmax(l, 1e-300)
elig[l < 1e-12] <- 1; elig <- pmax(elig, 1e-9)

windows <- c(0, 5, 5.5, 6, 8)
labels <- c("[0,5)", "[5,5.5)", "[5.5,6)", "[6,8)")
simulate_toy <- function(method) {
  parity <- integer(N); last <- rep(-Inf, N)
  nw <- length(windows) - 1
  births <- exposure <- lost <- matrix(0, nw, 3)
  for (i in 1:steps) {
    t <- tt[i]; w <- findInterval(t, windows)
    g <- group(parity); counts <- tabulate(g, 3)
    exposure[w, ] <- exposure[w, ] + counts * dt
    eligible <- (t - last) >= B - 1e-9
    m <- base(parity) * mult(t)
    h <- switch(method,
      none = m,
      composition = m / elig[cbind(i, pmin(parity, K) + 1)],
      runtime = {
        n_eligible <- tabulate(g[eligible], 3)
        lost[w, ] <- lost[w, ] + (n_eligible == 0) * counts * base(0:2) * mult(t) * dt
        m * (counts / pmax(n_eligible, 1))[g]
      })
    if (method == "none") eligible[] <- TRUE
    birth <- eligible & runif(N) < h * dt
    births[w, ] <- births[w, ] + tabulate(group(parity[birth]), 3)
    parity[birth] <- parity[birth] + 1L; last[birth] <- t
  }
  data.frame(method, window = rep(labels, 3), parity = rep(c("0", "1", "2+"), each = nw),
             births = c(births), exposure = c(exposure), rate = c(births / exposure),
             target = rep(base(0:2), each = nw) * rep(ifelse(windows[-1] <= 5, 1, 1.6), 3),
             expected_lost_births = c(lost))
}
results <- do.call(rbind, lapply(c("none", "composition", "runtime"), simulate_toy))
results$dt <- dt; results$women <- N
print(results, digits = 4, row.names = FALSE)
write.csv(results, "validation/milestone5/runtime-spacing-summary.csv", row.names = FALSE)
