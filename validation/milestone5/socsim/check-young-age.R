# Re-measure saved generational histories at age [18,19), years [80,100).
# Usage: Rscript check-young-age.R SCRATCH OUTPUT
args <- commandArgs(TRUE)
files <- list.files(args[1],pattern="persons.csv$",recursive=TRUE,full.names=TRUE)
files <- files[grepl("/long-",files)]
rows <- list()
for(file in files) {
  p <- read.csv(file)
  stopifnot(all(p$fmult[p$fem==1]==1))
  config <- read.csv(file.path(dirname(file),"rates.csv"))[1,]
  b <- split(p$dob[p$mom>0],p$mom[p$mom>0]); ix <- match(p$pid,names(b))
  exposure <- events <- numeric(5)
  for(i in which(p$fem==1 & p$dob+216<2401 & p$dob+228>2161)) {
    times <- if(is.na(ix[i])) numeric() else sort(b[[ix[i]]])
    lo <- max(2161,p$dob[i]+216)
    hi <- min(2401,p$dob[i]+228,if(p$dod[i]>0) p$dod[i]+1 else Inf)
    if(hi<=lo) next
    left <- c(p$dob[i],times+1);right <- c(times+1,Inf)
    parity <- pmin(seq_along(left)-1L,4L)+1L
    for(k in 1:5) exposure[k] <- exposure[k]+sum(pmax(0,
      pmin(right[parity==k],hi)-pmax(left[parity==k],lo)))/12
    selected <- which(times>=lo & times<hi)
    events <- events+tabulate(pmin(selected-1L,4L)+1L,nbins=5)
  }
  rows[[length(rows)+1L]] <- data.frame(bint=config$bint,seed=config$seed,
    time_start=80,time_end=100,age_start=18,age_end=19,
    parity=c("0","1","2","3","4+"),exposure,events,
    rate=ifelse(exposure>0,events/exposure,NA_real_),target=c(.2,.3,.4,.4,.4))
}
x <- do.call(rbind,rows)
stopifnot(nrow(x)==100L,!anyDuplicated(x[c("bint","seed","parity")]))
summary <- do.call(rbind,lapply(split(x,list(x$bint,x$parity)),function(z) {
  valid <- is.finite(z$rate)
  data.frame(bint=z$bint[1],parity=z$parity[1],n_valid=sum(valid),
    mean=mean(z$rate[valid]),mcse=sd(z$rate[valid])/sqrt(sum(valid)),target=z$target[1])
}))
write.csv(x,file.path(args[2],"age18-replications.csv"),row.names=FALSE)
write.csv(summary,file.path(args[2],"age18-summary.csv"),row.names=FALSE)
print(summary,row.names=FALSE)
