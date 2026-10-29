# Offline method-investigation diagnostic, not Hammel engine code or an
# acceptance test. Compares spacing corrections in an independent discrete-time
# toy simulation of an initialized cohort (all women parity zero at t = 0).
# Run from the repository root:
#   Rscript validation/milestone5/check-composition-spacing.R [dt]
# Writes validation/milestone5/composition-spacing-summary.csv (default dt only).
#
# Targets are conditional rates per woman-year at current parity (.2, .3, .4 for
# 2+), constant over ten years, spacing B = .75 years, no mortality.
# Corrections compared:
#   none        : no spacing (control) and raw targets
#   previous    : D099/SOCSIM-style h_k = m_k / (1 - B m_{k-1}), parity 0 unadjusted
#   composition : h_k = m_k / e_k(t), where e_k is the eligible fraction implied by
#                 the deterministic target cohort l_k(t):
#                 e_k = 1 - (inflow to parity k during the last B years) / l_k.
# See docs/parity-spacing-investigation.md for the derivation and limitations.
args <- commandArgs(trailingOnly = TRUE)
dt <- if (length(args)) as.numeric(args[[1]]) else 0.0025
set.seed(9601)
B <- 0.75; Tend <- 10; K <- 15; N <- 40000
m <- function(k) ifelse(k == 0, .2, ifelse(k == 1, .3, .4))
steps <- round(Tend / dt); tt <- (0:steps) * dt
mk <- m(0:K)

# Deterministic target cohort under the conditional (all-woman) target rates,
# tracking actual parity 0..K (K absorbing; negligible mass by t = 10).
l <- matrix(0, steps + 1, K + 1); l[1, 1] <- 1
for (i in 1:steps) {
  out <- mk * l[i, ]; out[K + 1] <- 0
  l[i + 1, ] <- l[i, ] - dt * out + dt * c(0, out[-(K + 1)])
}
inflow <- cbind(0, sweep(l[, -(K + 1)], 2, mk[-(K + 1)], `*`))
cumin <- apply(inflow, 2, cumsum) * dt
lagi <- pmax(0, (0:steps) - round(B / dt)) + 1
recent <- cumin - cumin[lagi, ]
elig <- 1 - recent / pmax(l, 1e-300)
elig[l < 1e-12] <- 1
elig <- pmax(elig, 1e-9)  # near-zero only where no simulated woman is eligible

simulate_toy <- function(method, spacing) {
  parity <- integer(N); last <- rep(-Inf, N)
  births <- exposure <- numeric(3)
  group <- function(p) pmin(p, 2L) + 1L
  for (i in 1:steps) {
    t <- tt[i]
    exposure <- exposure + tabulate(group(parity), 3) * dt
    h <- m(parity)
    if (method == "previous") h <- ifelse(parity == 0, h, h / (1 - spacing * m(parity - 1)))
    if (method == "composition") h <- h / elig[cbind(i, pmin(parity, K) + 1)]
    eligible <- (t - last) >= spacing - 1e-9
    birth <- eligible & runif(N) < h * dt
    births <- births + tabulate(group(parity[birth]), 3)
    parity[birth] <- parity[birth] + 1L; last[birth] <- t
  }
  data.frame(parity = c("0", "1", "2+"), births, exposure,
             rate = births / exposure, target = c(.2, .3, .4))
}

results <- rbind(
  cbind(correction = "none", spacing = 0, simulate_toy("none", 0)),
  cbind(correction = "previous", spacing = B, simulate_toy("previous", B)),
  cbind(correction = "composition", spacing = B, simulate_toy("composition", B)))
results$dt <- dt
results$women <- N
print(results, digits = 4, row.names = FALSE)
if (!length(args)) {
  write.csv(results, "validation/milestone5/composition-spacing-summary.csv",
            row.names = FALSE)
}
