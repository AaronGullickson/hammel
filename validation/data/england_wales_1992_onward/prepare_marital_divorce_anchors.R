# Offline derivation from transcribed ONS 2005 life table and retained 2010
# cumulative incidence. Run from the package root. No 1992 hazard is invented.
local({
  root <- "validation/data/england_wales_1992_onward"
  write_table <- function(x, name) {
    write.csv(tibble::as_tibble(x), file.path(root,name), row.names=FALSE,na="")
  }
  life <- read.csv(file.path(root,"marital_life_table_2005.csv"))
  stopifnot(nrow(life)==52L, identical(life$dur_start,c(0:50,55L)),
    all(life$dur_end>life$dur_start),
    all(life$prob_divorce>=0 & life$prob_death>=0),
    all(life$prob_divorce+life$prob_death<1),
    all(abs(head(life$surviving_marriages,-1)-head(life$divorce_events,-1)-
      head(life$death_events,-1)-tail(life$surviving_marriages,-1))<=1))
  # Columns g/h in the paper apply to the END of each interval. In particular
  # source duration 4 gives incidence by anniversary 5, not anniversary 4.
  # Integer life-table event counts retain more precision than the 0.1%-rounded
  # published cumulative columns (which hide small nonzero late divorces).
  cumulative <- data.frame(year=2005,anniversary=life$dur_end,
    cumulative_divorce=cumsum(life$divorce_events)/life$surviving_marriages[1],
    cumulative_death=cumsum(life$death_events)/life$surviving_marriages[1])
  stopifnot(all(abs(cumulative$cumulative_divorce-
    life$published_cumulative_divorce_end)<=.00051),
    all(abs(cumulative$cumulative_death-life$published_cumulative_death_end)<=.00051))
  write_table(cumulative,"marital_divorce_death_cumulative_2005.csv")

  q <- life[c("year","dur_start","dur_end","prob_divorce","prob_death")]
  write_table(q,"marital_divorce_death_probabilities_2005.csv")

  recent <- read.csv(file.path(root,"marital_divorce_death_cumulative_2010.csv"))
  anchor <- rbind(cumulative,recent)
  anchor <- anchor[anchor$anniversary %in% seq(5,50,5),]
  out <- do.call(rbind,lapply(split(anchor,anchor$year),function(x) {
    x <- x[order(x$anniversary),]
    stopifnot(length(x$anniversary)==10L,all(x$anniversary==seq(5,50,5)))
    previous_divorce <- c(0,head(x$cumulative_divorce,-1))
    previous_death <- c(0,head(x$cumulative_death,-1))
    surviving <- 1-previous_divorce-previous_death
    pd <- (x$cumulative_divorce-previous_divorce)/surviving
    pm <- (x$cumulative_death-previous_death)/surviving
    pe <- pd+pm
    stopifnot(all(pd>=0),all(pm>=0),all(pe<1))
    total <- -log1p(-pe)/5
    divorce <- ifelse(pe>0,total*pd/pe,0)
    mortality <- ifelse(pe>0,total*pm/pe,0)
    # Validate both causes against the original conditional interval outcomes.
    recovered <- -expm1(-5*(divorce+mortality))
    stopifnot(max(abs(ifelse(pe>0,recovered*divorce/total,0)-pd))<1e-12,
      max(abs(ifelse(pe>0,recovered*mortality/total,0)-pm))<1e-12)
    data.frame(year=x$year,union_type="marriage",dur_start=x$anniversary-5,
      dur_end=ifelse(x$anniversary==50,Inf,x$anniversary),
      observed_dur_end=x$anniversary,rate=divorce)
  }))
  rownames(out) <- NULL
  write_table(out,"marital_divorce_hazard_anchors_2005_2010.csv")
  message("Rebuilt 2005/2010 anchors; 1993-94 remains benchmarks only, not a complete hazard schedule.")
})
