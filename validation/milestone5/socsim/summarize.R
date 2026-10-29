# Run after all fresh-process jobs: Rscript summarize.R SCRATCH OUTPUT
args <- commandArgs(TRUE)
files <- list.files(args[1],pattern="rates.csv$",recursive=TRUE,full.names=TRUE)
stopifnot(length(files)==40L)
x <- do.call(rbind,lapply(files,read.csv))
x$experiment <- ifelse(x$months==120,"initialized_cohort","generational")
stopifnot(nrow(x)==360L, !anyDuplicated(x[c("experiment","bint","seed","time_start","parity")]))
summary <- do.call(rbind,lapply(split(x,list(x$experiment,x$bint,x$time_start,x$parity),drop=TRUE),function(z) {
  valid <- is.finite(z$rate)
  data.frame(experiment=z$experiment[1],bint=z$bint[1],time_start=z$time_start[1],
    time_end=z$time_end[1],parity=z$parity[1],n_valid=sum(valid),
    mean=mean(z$rate[valid]),mcse=sd(z$rate[valid])/sqrt(sum(valid)),
    target=z$target[1])
}))
write.csv(x,file.path(args[2],"replications.csv"),row.names=FALSE)
write.csv(summary,file.path(args[2],"summary.csv"),row.names=FALSE)
print(subset(summary,experiment=="initialized_cohort" | time_start==80),row.names=FALSE)
