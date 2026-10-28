# Offline D086 parameter recovery with different marriage/cohabitation models.
# Rscript validation/validate-cohabitation-matching.R LIB_OR_dev PREFIX [REPS]
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=2L)
if (args[1]=="dev") pkgload::load_all(compile=FALSE,quiet=TRUE) else library(hammel,lib.loc=args[1])
suppressPackageStartupMessages(library(survival))
reps <- if(length(args)>2L) as.integer(args[3]) else 10L
choice_data <- function(result,type,n=1500L) {
  p <- result$persons; u <- result$unions
  wife <- match(u$wife_id,p$person_id); husband <- match(u$husband_id,p$person_id)
  male <- which(p$sex=="m")
  eligible_rows <- which(u$union_type==type & is.na(u$previous_cohab_id))
  selected <- eligible_rows[sample.int(length(eligible_rows),min(n,length(eligible_rows)))]
  do.call(rbind,lapply(seq_along(selected),function(k) {
    row <- selected[k]; now <- u$start_time[row]; chosen <- husband[row]
    partnered <- husband[u$start_time<=now & (is.na(u$end_time)|u$end_time>now)]
    partnered <- partnered[partnered!=chosen]
    pool <- male[p$birth_time[male]+15<=now &
      (is.na(p$death_time[male])|p$death_time[male]>now) & !male %in% partnered]
    stopifnot(chosen %in% pool)
    other <- pool[pool!=chosen]
    candidates <- c(chosen,other[sample.int(length(other),min(50L,length(other)))])
    gf <- p$group[wife[row]]; gm <- p$group[candidates]
    data.frame(set=k,chosen=seq_along(candidates)==1L,
      age_diff=p$birth_time[wife[row]]-p$birth_time[candidates],
      AB=as.numeric(gf=="A" & gm=="B"),BA=as.numeric(gf=="B" & gm=="A"),
      log_weight=ifelse(gm=="A",u$.log_weight_A[row],u$.log_weight_B[row]))
  }))
}
rows <- list()
for (replication in seq_len(reps)) {
  set.seed(286000L+replication)
  pop <- create_baseline_population(20000,50,2000,group_weights=c(A=3,B=1))
  search <- data.frame(age_start=c(15,20,25,30,35,45,60),age_end=c(20,25,30,35,45,60,Inf),
    rate=c(.02,.08,.14,.1,.06,.03,.005))
  search <- rbind(transform(search,previous_union="none"),
    transform(search,previous_union="dissolved",rate=rate*1.4),
    transform(search,previous_union="widowed",rate=rate*.7))
  prob <- expand.grid(group=c("A","B"),previous_union=c("none","dissolved","widowed"),
    age_start=c(15,25,35,50),stringsAsFactors=FALSE)
  prob$age_end <- c(25,35,50,Inf)[match(prob$age_start,c(15,25,35,50))]
  prob$prob <- c(.7,.6,.4,.25)[match(prob$age_start,c(15,25,35,50))] *
    ifelse(prob$group=="B",.9,1) * ifelse(prob$previous_union=="none",1,.8)
  model <- hammel_model() |> set_population(pop) |>
    set_mortality(data.frame(sex=c("m","f"),age_start=0,age_end=Inf,rate=.01)) |>
    set_fertility(data.frame(age_start=15,age_end=50,
      union_type=c("unpartnered","cohab","marriage"),rate=c(.04,.10,.12))) |>
    set_union_formation("rate",rates=search,
      matching=create_age_gap_matching_model(2,5,group_exogamy=c("A:B"=-3)),
      cohab_matching=create_age_gap_matching_model(0,7,group_exogamy=c("A:B"=-1)),
      cohab_prob=prob,union_convert_rate=data.frame(dur_start=c(0,1,5),dur_end=c(1,5,Inf),rate=c(0,.3,.15)),
      group_prop_hazard=c(B=1.2),diagnostics=TRUE,kinship_generations=0L) |>
    set_union_dissolution(data.frame(union_type=c("cohab","marriage"),dur_start=0,dur_end=Inf,rate=c(.08,.02)))
  result <- simulate(model,60,287000L+replication)
  for (type in c("marriage","cohab")) {
    set.seed(288000L+replication)
    data <- choice_data(result,type)
    fit <- clogit(chosen~age_diff+I(age_diff^2)+AB+BA+offset(log_weight)+strata(set),data=data,method="efron")
    beta <- coef(fit)
    rows[[length(rows)+1L]] <- data.frame(replication,union_type=type,
      choices=length(unique(data$set)),AB=beta[["AB"]],BA=beta[["BA"]],
      log_odds=beta[["AB"]]+beta[["BA"]],linear=beta[["age_diff"]],quadratic=beta[["I(age_diff^2)"]],
      target=if(type=="marriage") -6 else -2,conversions=result$simulation$conversions)
    message(type," rep=",replication," log OR=",round(tail(rows,1)[[1]]$log_odds,4))
  }
}
results <- do.call(rbind,rows)
summary <- do.call(rbind,lapply(split(results,results$union_type),function(x)
  do.call(rbind,lapply(c("log_odds","AB","BA","linear","quadratic"),function(parameter)
    data.frame(union_type=x$union_type[1],parameter,mean=mean(x[[parameter]]),mcse=sd(x[[parameter]])/sqrt(nrow(x)))))))
write.csv(results,paste0(args[2],"-replications.csv"),row.names=FALSE)
write.csv(summary,paste0(args[2],"-summary.csv"),row.names=FALSE)
saveRDS(list(results=results,summary=summary,session=sessionInfo(),source_md5=tools::md5sum("src/simulate.cpp")),paste0(args[2],".rds"))
print(summary,row.names=FALSE)
