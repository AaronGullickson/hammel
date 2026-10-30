# Offline diagnostic for D099's SOCSIM-style approximation, not an acceptance test.
# Run from the repository root: Rscript validation/milestone5/check-parity-spacing.R
# All founders are women aged 20 at parity zero; no mortality or birth descendants
# within ten years. These are finite initialized cohorts, not stable populations.
pkgload::load_all(quiet=TRUE)
results <- list()
for (spacing in c(0,.75)) for (replication in 1:10) {
  n <- 3000L
  people <- data.frame(person_id=seq_len(n),sex="f",birth_time=-20,
    death_time=NA_real_,mother_id=NA_integer_,father_id=NA_integer_)
  rates <- data.frame(age_start=15,age_end=Inf,parity_start=0:2,
    parity_end=c(1,2,Inf),rate=c(.2,.3,.4))
  model <- hammel::hammel_model() |>
    hammel::set_population(hammel::hammel_population(people,end_time=0)) |>
    hammel::set_mortality(data.frame(sex=c("m","f"),age_start=0,age_end=Inf,rate=0)) |>
    hammel::set_fertility(rates,birth_interval=spacing) |>
    hammel::set_union_formation("on_demand",cancel_birth_if_no_match=FALSE)
  population <- hammel::simulate(model,10,9500+replication)
  born <- population$persons[population$persons$mother_id %in% seq_len(n), ]
  times <- split(born$birth_time,born$mother_id)
  exposure <- events <- numeric(3)
  for (i in seq_len(n)) {
    births <- sort(times[[as.character(i)]])
    intervals <- diff(c(0,births,10))
    parity <- pmin(seq_along(intervals)-1L,2L)+1L
    for (k in 1:3) exposure[k] <- exposure[k]+sum(intervals[parity==k])
    events <- events+tabulate(pmin(seq_along(births)-1L,2L)+1L,nbins=3)
  }
  results[[length(results)+1L]] <- data.frame(spacing,replication,
    parity=c("0","1","2+"),exposure,events,rate=events/exposure)
}
results <- do.call(rbind,results)
summary <- do.call(rbind,lapply(split(results,list(results$spacing,results$parity)),function(x) {
  data.frame(spacing=x$spacing[1],parity=x$parity[1],mean=mean(x$rate),
    mcse=sd(x$rate)/sqrt(nrow(x)),target=c('0'=.2,'1'=.3,'2+'=.4)[x$parity[1]])
}))
write.csv(results,"validation/milestone5/parity-spacing-replications.csv",row.names=FALSE)
write.csv(summary,"validation/milestone5/parity-spacing-summary.csv",row.names=FALSE)
print(summary,row.names=FALSE)
