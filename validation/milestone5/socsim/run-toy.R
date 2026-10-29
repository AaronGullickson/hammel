# Independent SOCSIM reference experiment; no engine source copied or modified.
# Usage: Rscript run-toy.R LIBRARY SCRATCH BINT SEED [MONTHS] [WOMEN]
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
library(rsocsim)
folder <- normalizePath(args[2], mustWork=FALSE)
dir.create(folder, recursive=TRUE, showWarnings=FALSE)
bint <- as.integer(args[3]); seed <- as.integer(args[4])
months <- if(length(args)>4) as.integer(args[5]) else 120L
n <- if(length(args)>5) as.integer(args[6]) else 3000L
# Last observed date 1200 anchors simulation months 1201 onward.
p <- data.frame(pid=seq_len(n+1L), fem=c(rep(1L,n),0L), group=1,
  nev=0,dob=c(rep(960,n),1200),mom=0,pop=0,nesibm=0,nesibp=0,
  lborn=0,marid=0,mstat=1,dod=0,fmult=1)
write.table(p,file.path(folder,"initial.opop"),row.names=FALSE,col.names=FALSE)
file.create(file.path(folder,"initial.omar"))
# Fertility inputs are monthly hazard parameters (load.cpp::add_rate_block).
# Founders have multiplier one; alpha=1 and betaT=1 also fix all descendants
# at one, independently of the heterogeneity toggle.
# Exactly zero input is internally floored by SOCSIM, negligibly here.
lines <- c("hetfert 0","alpha 1","betaT 1","random_father 0","marriage_after_childbirth 0",
  "marriage_queues 1",paste("bint",bint),"segments 1","input_file initial",
  paste("duration",months))
for(parity in 0:2) lines <- c(lines,paste("birth 1 F single",parity),
  "15 0 0",paste("40 0",format(c(.2,.3,.4)[parity+1L]/12,digits=17)),"120 0 0")
for(sex in c("F","M")) lines <- c(lines,paste("death 1",sex,"single"),
  "60 0 0","120 0 0.1",paste("marriage 1",sex,"single"),"120 0 0")
writeLines(c(lines,"run"),file.path(folder,"toy.sup"))
status <- socsim(folder,"toy.sup",seed=seed,suffix="toy")
stopifnot(status==1)
out <- read_opop(folder=folder,supfile="toy.sup",seed=seed,suffix="toy",quiet=TRUE)
stopifnot(all(out$fmult[out$fem==1]==1), all(out$mstat==1))
write.csv(out,file.path(folder,"persons.csv"),row.names=FALSE)
# Independent histories and monthly exposure, not SOCSIM's own rate summaries.
# A woman's event-month contributes pre-birth parity exposure; new parity starts
# the following month. This matches one monthly person-time risk census.
births <- split(out$dob[out$mom>0],out$mom[out$mom>0])
birth_index <- match(out$pid,names(births))
windows <- if(months==120) c(0,120) else seq(0,months,by=240)
rows <- list()
for(w in seq_len(length(windows)-1L)) {
  exposure <- events <- numeric(3)
  for(i in which(out$fem==1 & (months>120 | out$pid<=n))) {
    b <- if(is.na(birth_index[i])) numeric() else sort(births[[birth_index[i]]])
    lo <- max(1201+windows[w],out$dob[i]+if(months>120) 20*12 else 0)
    hi <- min(1201+windows[w+1L],if(out$dod[i]>0) out$dod[i]+1 else Inf,
      if(months>120) out$dob[i]+35*12 else Inf)
    if(hi<=lo) next
    left <- c(out$dob[i],b+1); right <- c(b+1,Inf)
    parity <- pmin(seq_along(left)-1L,2L)+1L
    for(k in 1:3) exposure[k] <- exposure[k]+sum(pmax(0,
      pmin(right[parity==k],hi)-pmax(left[parity==k],lo)))/12
    selected <- which(b>=lo & b<hi)
    events <- events+tabulate(pmin(selected-1L,2L)+1L,nbins=3)
  }
  rows[[w]] <- data.frame(bint,seed,months,n,time_start=windows[w]/12,
    time_end=windows[w+1L]/12,parity=c("0","1","2+"),exposure,events,
    rate=ifelse(exposure>0,events/exposure,NA_real_),target=c(.2,.3,.4))
}
result <- do.call(rbind,rows)
if(months==120) stopifnot(abs(sum(result$exposure)-n*10)<1e-7,
  sum(result$events)==sum(out$mom %in% seq_len(n)))
write.csv(result,file.path(folder,"rates.csv"),row.names=FALSE)
print(result)
