# Analytical fixed-market calculations; does not load or change Hammel.
# Run: Rscript validation/matching_column_adjustment.R
# Output: validation/matching_column_adjustment/ (or HAMMEL_COLUMN_OUTPUT).

log_sum <- function(x) { z <- max(x); z + log(sum(exp(x - z))) }
normalize <- function(x) exp(x - log_sum(x))

# Matrices are internal calculation objects; output tables are tibbles.
calculate_market <- function(f, m, utility, alpha, additive = FALSE) {
  stopifnot(length(f) == nrow(utility), length(m) == ncol(utility),
    all(f > 0), all(m > 0), abs(sum(f) - 1) < 1e-12,
    abs(sum(m) - 1) < 1e-12, all(is.finite(utility)), is.finite(alpha))
  log_columns <- vapply(seq_along(m), function(j)
    log_sum(log(f) + utility[, j]), numeric(1))
  prob <- matrix(0, length(f), length(m))
  for (i in seq_along(f)) {
    if (additive) {
      # Equal-sized sampled pools; summed individual male-to-ego probabilities.
      a <- normalize(log(m) + utility[i, ])
      b <- exp(log(m) + utility[i, ] - log_columns)
      prob[i, ] <- (a + b) / sum(a + b)
    } else prob[i, ] <- normalize(log(m) + utility[i, ] - alpha * log_columns)
  }
  table <- prob * f
  stopifnot(max(abs(rowSums(prob) - 1)) < 1e-12,
    max(abs(rowSums(table) - f)) < 1e-12, abs(sum(table) - 1) < 1e-12)
  list(table = table, probability = prob, log_columns = log_columns)
}

summarize_two <- function(table) {
  tibble::tibble(AA = table[1, 1], AB = table[1, 2],
    BA = table[2, 1], BB = table[2, 2],
    log_odds = log(table[1, 2]) + log(table[2, 1]) -
      log(table[1, 1]) - log(table[2, 2]),
    gender_ratio = table[1, 2] / table[2, 1])
}
make_utility <- function(ab, ba = ab) matrix(c(0, ba, ab, 0), 2)
fit_two_alpha <- function(p, beta) {
  if (abs(p - .5) < 1e-12 || abs(beta) < 1e-12) return(NA_real_)
  objective <- function(alpha) {
    t <- calculate_market(c(p, 1-p), c(p, 1-p), make_utility(beta), alpha)$table
    log(t[1, 2]) - log(t[2, 1])
  }
  uniroot(objective, c(-32, 32), tol = 1e-12)$root
}
# Independent closed-form two-group solution, with t = b_B / b_A.
closed_two_alpha <- function(p, beta) {
  if (abs(p - .5) < 1e-12 || abs(beta) < 1e-12) return(NA_real_)
  q <- exp(beta)
  t <- (-(2*p-1)*q + sqrt((2*p-1)^2*q^2 + 4*p*(1-p))) / (2*(1-p))
  log(t) / log((p + (1-p)*q) / ((1-p) + p*q))
}

reference_alpha <- fit_two_alpha(.75, -3)
rows <- list()
for (p in c(.01, .05, .1, .25, .5, .75, .9, .95, .99)) {
  for (beta in c(-6, -3, -1, -.25, 0, .25, 1, 3)) {
    fitted <- fit_two_alpha(p, beta)
    closed <- closed_two_alpha(p, beta)
    stopifnot(is.na(fitted) || abs(fitted - closed) < 1e-9)
    methods <- c(female_choice = 0, geometric = .5,
      fixed_reference = reference_alpha, case_fitted = if (is.na(fitted)) .5 else fitted)
    for (method in names(methods)) {
      alpha <- methods[[method]]
      table <- calculate_market(c(p, 1-p), c(p, 1-p), make_utility(beta), alpha)$table
      result <- summarize_two(table)
      stopifnot(abs(result$log_odds - 2*beta) < 1e-10)
      if (method == 'case_fitted') stopifnot(abs(log(result$gender_ratio)) < 1e-9)
      rows[[length(rows)+1L]] <- cbind(tibble::tibble(p_A = p, beta,
        method, alpha, fitted_alpha = fitted), result)
    }
  }
}
two <- tibble::as_tibble(do.call(rbind, rows))

# Expanded 100-by-100 individual market checks group aggregation and weighting.
group <- c(rep(1L, 75), rep(2L, 25))
for (alpha in c(0, .5, reference_alpha)) {
  u <- make_utility(-3)
  expanded <- calculate_market(rep(.01, 100), rep(.01, 100), u[group, group], alpha)$table
  collapsed <- matrix(0, 2, 2)
  for (i in 1:2) for (j in 1:2)
    collapsed[i, j] <- sum(expanded[group == i, group == j])
  expected <- calculate_market(c(.75, .25), c(.75, .25), u, alpha)$table
  stopifnot(max(abs(collapsed - expected)) < 1e-12)
}
additive <- summarize_two(calculate_market(c(.75,.25), c(.75,.25),
  make_utility(-3), 0, additive=TRUE)$table)
stopifnot(abs(additive$log_odds - (-5.770809024484435)) < 1e-12)

# Freeze the one reference exponent; do not fit equality to asymmetric cases.
cases <- tibble::tibble(scenario = c('symmetric_reference', 'more_B_women',
  'more_B_men', 'directed_AB_easier', 'directed_BA_easier'),
  female_A = c(.75,.65,.85,.75,.75), male_A = c(.75,.85,.65,.75,.75),
  beta_AB = c(-3,-3,-3,-1,-5), beta_BA = c(-3,-3,-3,-5,-1))
rows <- list()
for (k in seq_len(nrow(cases))) for (method in c('female_choice','geometric','fixed_reference')) {
  case <- cases[k, ]; alpha <- c(female_choice=0, geometric=.5, fixed_reference=reference_alpha)[[method]]
  t <- calculate_market(c(case$female_A,1-case$female_A),
    c(case$male_A,1-case$male_A), make_utility(case$beta_AB,case$beta_BA),alpha)$table
  result <- summarize_two(t)
  stopifnot(abs(result$log_odds - case$beta_AB - case$beta_BA) < 1e-10)
  rows[[length(rows)+1L]] <- cbind(case,tibble::tibble(method,alpha),result)
}
asymmetric <- tibble::as_tibble(do.call(rbind,rows))

# Three symmetric groups: a scalar need not balance all flows simultaneously.
# Fit each pair separately, then find the scalar minimizing all log imbalances.
rows <- list()
for (shares in list(c(.6,.3,.1),c(.8,.15,.05),rep(1/3,3))) {
  for (scenario in c('common_effect','different_pair_effects')) {
    u <- matrix(-3,3,3); diag(u) <- 0
    if (scenario == 'different_pair_effects') {
      u[1,3] <- u[3,1] <- -1; u[2,3] <- u[3,2] <- -2
    }
    imbalance <- function(alpha) {
      t <- calculate_market(shares,shares,u,alpha)$table
      c(log(t[1,2]/t[2,1]),log(t[1,3]/t[3,1]),log(t[2,3]/t[3,2]))
    }
    best <- optimize(function(alpha)sum(imbalance(alpha)^2),c(-8,8),tol=1e-12)$minimum
    for (method in c('geometric','fixed_reference','best_scalar')) {
      alpha <- c(geometric=.5,fixed_reference=reference_alpha,best_scalar=best)[[method]]
      t <- calculate_market(shares,shares,u,alpha)$table
      for (pair in list(c(1,2),c(1,3),c(2,3))) {
        i <- pair[1];j <- pair[2]
        odds <- log(t[i,j])+log(t[j,i])-log(t[i,i])-log(t[j,j])
        stopifnot(abs(odds-u[i,j]-u[j,i])<1e-10)
        rows[[length(rows)+1L]] <- tibble::tibble(shares=paste(shares,collapse='/'),
          scenario,method,alpha,pair=paste(LETTERS[pair],collapse=''),
          gender_ratio=t[i,j]/t[j,i],log_odds=odds,target=u[i,j]+u[j,i])
      }
    }
  }
}
three <- tibble::as_tibble(do.call(rbind,rows))
output <- Sys.getenv('HAMMEL_COLUMN_OUTPUT','validation/matching_column_adjustment')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
for (name in c('two','asymmetric','three','additive'))
 write.csv(get(name),file.path(output,paste0(name,'.csv')),row.names=FALSE)
saveRDS(list(time=Sys.time(),session=sessionInfo(),reference_alpha=reference_alpha,
 script_hash=tools::md5sum('validation/matching_column_adjustment.R')),
 file.path(output,'metadata.rds'))
cat('Reference alpha:',format(reference_alpha,digits=12),'\n')
print(two[two$method=='case_fitted' & two$p_A %in% c(.5,.75,.95) &
 two$beta %in% c(-3,-1,.25),c('p_A','beta','fitted_alpha','log_odds','gender_ratio')],n=30)
print(asymmetric[asymmetric$method=='fixed_reference',c('scenario','log_odds','gender_ratio')])
print(three[three$shares=='0.6/0.3/0.1' & three$method=='best_scalar',],n=30)
cat('All analytical invariants passed. Outputs:',output,'\n')
