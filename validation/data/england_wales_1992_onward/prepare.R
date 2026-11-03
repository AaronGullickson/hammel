# Run from the package root. Offline extraction of published observations, not
# a calibrated simulation scenario. Requires readxl and tibble already installed.
root <- "validation/data/england_wales_1992_onward"
read_sheet <- function(file, sheet) {
  as.data.frame(readxl::read_excel(file.path(root, file), sheet = sheet,
    col_names = FALSE, .name_repair = "minimal"), check.names = FALSE)
}
number <- function(x) suppressWarnings(as.numeric(x))
save_table <- function(x, name) {
  write.csv(tibble::as_tibble(x), file.path(root, name), row.names = FALSE, na = "")
}
age_bounds <- function(label) {
  nums <- number(regmatches(label, gregexpr("[0-9]+", label))[[1]])
  if (length(nums) == 2L) return(c(nums[1], nums[2] + 1))
  if (length(nums) == 1L && grepl("over|\\+", label)) return(c(nums, Inf))
  stop("Unrecognised age label: ", label)
}

# Published Table B1. The end bands describe denominator ages, not a claim
# that the birth numerators exclude births below 16 or above 44. Preserve both.
fertility <- rbind(
  cbind(year = 1992, rbind(c(286.7,192,22), c(196.3,127.9,31),
    c(169.5,107.9,32.7), c(103.3,101.2,26.3), c(35.5,62.3,13.3), c(5.7,16.6,2.9))),
  cbind(year = 2001, rbind(c(287.1,269.8,21.1), c(228.9,128.2,26.5),
    c(164.9,90.4,23.6), c(117,95.3,19), c(46.6,65.5,11.5), c(8.6,19.8,3.6))),
  cbind(year = 2007, rbind(c(456.3,330.5,18.6), c(286.5,144,28.7),
    c(219.3,108.8,28.4), c(164.5,114.4,24.1), c(68.3,89.1,15.6), c(12.3,29.2,4.8))))
f <- do.call(rbind, lapply(1:3, function(j) data.frame(
  year = fertility[,1], source_age_label = rep(c("Under 20","20-24","25-29",
    "30-34","35-39","40 and over"), 3),
  denominator_age_start = rep(c(16,20,25,30,35,40),3),
  denominator_age_end = rep(c(20,25,30,35,40,45),3),
  union_type = c("marriage","cohab","unpartnered")[j],
  births_per_1000_women = fertility[,j+1], rate = fertility[,j+1]/1000)))
stopifnot(nrow(f) == 54L, all(f$rate >= 0))
save_table(f, "fertility_by_union_type_1992_2001_2007.csv")
save_table(subset(f, year == 1992), "fertility_by_union_type_1992.csv")

# Table A4: source estimates rounded to hundreds, so sums need not equal totals.
stock <- rbind(
  cbind(year=1992,rbind(c(27200,79300,1128600),c(439200,328900,1128600),
    c(1114800,322100,645200),c(1334800,184600,390100),c(1289700,112600,293200),
    c(1379300,90500,285900))),
  cbind(year=2001,rbind(c(16200,60100,1104300),c(178000,318700,1027200),
    c(625400,440400,721200),c(1142400,349500,626500),c(1356300,256400,565600),
    c(1292900,175100,466500))),
  cbind(year=2007,rbind(c(6800,51500,1327500),c(135600,388300,1256200),
    c(486200,496800,772600),c(839600,347700,547400),c(1191300,284800,551600),
    c(1334100,216200,554200))))
s <- f[c("year","source_age_label","denominator_age_start","denominator_age_end","union_type")]
s$women <- unlist(lapply(2:4,function(j) stock[,j]),use.names=FALSE)
save_table(s, "female_union_type_population_1992_2001_2007.csv")

# Annual female population denominators: old demographic estimates through
# 2001; survey-based revised estimates 2002 onward. Keep the method break.
pop <- list()
for (year in 1992:2001) {
  d <- read_sheet("ons_marital_population_1971_2001.xls", as.character(year))
  rows <- which(grepl("^[0-9]+$|^[0-9]+ and over$", d[[1]]))
  for (j in seq_along(rows)) {
    i <- rows[j]
    a <- number(sub(" and over", "", d[i,1]))
    for (status in c("Single","Widowed","Divorced")) {
      col <- c(Single=9, Widowed=11, Divorced=12)[status]
      pop[[length(pop)+1L]] <- data.frame(year=year, marital_status=status,
        age_start=a, age_end=if (grepl("over",d[i,1])) Inf else a+1,
        population=number(d[i,col])*1000, cv_code=NA_character_,
        population_ci_half_width=NA_real_, denominator_method="historical_demographic_estimate")
    }
  }
}
d <- read_sheet("ons_marital_population_2002_2019.xlsx", "1 E&W Marital Status")
start <- which(d[[1]] == "FEMALES")
hdr <- start+1L
year_cols <- which(number(unlist(d[hdr,])) %in% 2002:2019)
for (i in seq.int(hdr+1L,nrow(d))) {
  status <- d[i,1]
  label <- d[i,2]
  if (is.na(status) || !status %in% c("Single","Widowed","Divorced") ||
      is.na(label) || !grepl("^[0-9]+[-+]|^[0-9]+$",label)) next
  bounds <- age_bounds(label)
  for (j in year_cols) {
    pop[[length(pop)+1L]] <- data.frame(year=number(d[hdr,j]),marital_status=status,
      age_start=bounds[1],age_end=bounds[2],population=number(d[i,j]),
      cv_code=as.character(d[i,j+1]),population_ci_half_width=number(d[i,j+2]),
      denominator_method="survey_based_2020_release")
  }
}
pop <- do.call(rbind,pop)
stopifnot(setequal(unique(pop$year),1992:2019), !anyDuplicated(pop[1:4]))
save_table(pop,"female_marital_population_1992_2019.csv")

# Counts and denominators use exactly the same age bands. No interpolation,
# population smoothing, imputation of suppressed cells, or split of open ages.
d <- read_sheet("ons_marriages_2019.xlsx","3a")
header <- which(d[[1]] == "Year")
cols <- which(grepl("^Women",unlist(d[header,])) &
  grepl("Under 20|[0-9]+ to [0-9]+|80 and over",unlist(d[header,])))
out <- list()
for (year in 1992:2019) {
  for (j in cols) {
    label <- sub("Women[\r\n]+","",d[header,j])
    bounds <- if(label == "Under 20") c(16,20) else age_bounds(label)
    for (kind in c("first_marriage","remarriage","all_marriages")) {
      statuses <- switch(kind,first_marriage="Single",remarriage=c("Divorced","Widowed"),
        all_marriages=c("Single","Divorced","Widowed"))
      rows <- which(d[[1]] == as.character(year) & d[[2]] %in% statuses)
      count <- sum(number(unlist(d[rows,j])))
      p <- pop[pop$year == year & pop$marital_status %in% statuses &
        pop$age_start >= bounds[1] & pop$age_end <= bounds[2], , drop=FALSE]
      stopifnot(nrow(p) > 0L)
      denom <- sum(p$population) # NA stays NA; a suppressed cell is not zero.
      out[[length(out)+1L]] <- data.frame(year=year,age_start=bounds[1],age_end=bounds[2],
        marriage_order=kind,marriages=count,female_midyear_population=denom,
        rate=if(is.finite(denom) && denom>0) count/denom else NA_real_,
        rate_status=if(!is.finite(denom)) "missing_population_cell" else if(denom==0)
          "zero_population_denominator" else "derived_central_rate",
        denominator_method=unique(p$denominator_method))
    }
  }
}
out <- do.call(rbind,out)
stopifnot(nrow(out)==28L*12L*3L, all(out$rate >= 0,na.rm=TRUE))
# Counts by band must reconcile with published status-specific age-known totals.
for (year in 1992:2019) {
  for (kind in c("first_marriage","remarriage","all_marriages")) {
    statuses <- switch(kind,first_marriage="Single",remarriage=c("Divorced","Widowed"),
      all_marriages=c("Single","Divorced","Widowed"))
    rows <- which(d[[1]]==as.character(year) & d[[2]] %in% statuses)
    unstated <- number(d[rows,40])
    # ONS note 10: age-not-stated cases were distributed over ages in 1961-93.
    if (year <= 1993) unstated[] <- 0
    expected <- sum(number(d[rows,22]))-sum(unstated)-sum(number(d[rows,23]))
    stopifnot(sum(out$marriages[out$year==year & out$marriage_order==kind])==expected)
  }
}
save_table(out,"female_marriage_rates_1992_2019.csv")

# Also retain the published pooled rates; these have different native bands
# and remain useful where small status-specific survey cells are suppressed.
d <- read_sheet("ons_marriages_2019.xlsx","5a")
header <- which(d[[1]]=="Year")
cols <- which(grepl("^Women",unlist(d[header,])) &
  !grepl("Ages 16 and over",unlist(d[header,])))
rows <- which(number(d[[1]]) %in% 1992:2019)
out <- list()
for(i in rows) for(j in cols) {
  label <- sub("Women[\r\n]+","",d[header,j])
  bounds <- if(label=="Under 20") c(16,20) else age_bounds(label)
  out[[length(out)+1L]] <- data.frame(year=number(d[i,1]),
    age_start=bounds[1],age_end=bounds[2],rate=number(d[i,j])/1000)
}
out <- do.call(rbind,out)
stopifnot(nrow(out)==28L*11L, all(is.finite(out$rate)), all(out$rate>=0))
save_table(out,"female_published_pooled_marriage_rates_1992_2019.csv")

# Living arrangements are population stocks, never cohabitation-entry probabilities.
d <- read_sheet("ons_marital_population_2002_2019.xlsx","2 E&W Living Arrangements")
start <- which(d[[1]]=="FEMALES")
hdr <- start+1L
year_cols <- which(number(unlist(d[hdr,])) %in% 2002:2019)
out <- list()
for (i in seq.int(hdr+1L,nrow(d))) {
  label <- d[i,2]
  if(is.na(label) || !grepl("^[0-9]+[-+]|^[0-9]+$",label)) next
  bounds <- age_bounds(label)
  for(j in year_cols) out[[length(out)+1L]] <- data.frame(year=number(d[hdr,j]),
    age_start=bounds[1],age_end=bounds[2],living_arrangement=d[i,1],
    women=number(d[i,j]),cv_code=d[i,j+1],population_ci_half_width=number(d[i,j+2]))
}
save_table(do.call(rbind,out),"female_living_arrangements_2002_2019.csv")

# ONS Table 11a uses shared addresses at marriage, conditional on marrying.
# Preserve the all-marriages measure, including its source age labels.
d <- read_sheet("ons_marriages_2019.xlsx","11a")
rows <- which(number(d[[1]]) %in% 1994:2019)
save_table(data.frame(year=number(d[rows,1]),source_age_label=d[rows,2],
  proportion_cohabiting_before_marriage=number(d[rows,6])/100),
  "female_premarital_cohabitation_1994_2019.csv")

d <- read_sheet("ons_divorce_probability_2010.xls","Sheet")
rows <- which(is.finite(number(d[[1]])))
q <- data.frame(year=2010,dur_start=number(d[rows,1]),
  dur_end=number(d[rows,1])+1,prob_divorce=number(d[rows,2])/100)
stopifnot(all(q$prob_divorce>=0 & q$prob_divorce<1), !anyDuplicated(q$dur_start))
save_table(q,"marital_divorce_probabilities_by_duration_2010.csv")
d <- read.csv(file.path(root,"ons_divorce_death_cumulative_2010.csv"),skip=6)
save_table(data.frame(year=2010,anniversary=number(d[[1]]),
  cumulative_divorce=number(d[[2]])/100,cumulative_death=number(d[[3]])/100),
  "marital_divorce_death_cumulative_2010.csv")

# Cohort cumulative incidence, not period duration hazards. Retain 1992 onward.
d <- read_sheet("ons_divorces_2019.xlsx","Table 6")
yr <- number(d[[1]])
rows <- which(is.finite(yr) & yr>=1992 & yr<=2018)
hdr <- which(number(d[[3]]) == 1 & number(d[[12]]) == 10)[1]
anniv <- number(unlist(d[hdr,]))
cols <- which(is.finite(anniv) & anniv>=1 & anniv<=60)
out <- do.call(rbind,lapply(rows,function(i) data.frame(marriage_year=yr[i],
  marriages=number(d[i,2]),anniversary=anniv[cols],
  cumulative_divorce=number(unlist(d[i,cols]))/100)))
out <- subset(out,is.finite(cumulative_divorce))
stopifnot(nrow(out)>0, all(out$cumulative_divorce>=0 & out$cumulative_divorce<=1))
save_table(out,"marriage_cohort_cumulative_divorce_1992_2018.csv")

# Beaujouan & Ni Bhrolchain Table 3: women <45 at cohabitation start,
# competing-risk cumulative probabilities. Earliest cohort spans 1992.
c <- data.frame(cohort_start=c(1990,1995,2000),cohort_end=c(1995,2000,2005),
  marriage_5=c(.463,.399,.395),separation_5=c(.232,.305,.354),
  continuing_5=c(.305,.296,.252),marriage_10=c(.561,.494,NA),
  separation_10=c(.320,.381,NA),continuing_10=c(.120,.124,NA))
stopifnot(all(abs(rowSums(c[3:5])-1)<.002))
save_table(c,"great_britain_cohabitation_cohort_outcomes.csv")
p <- data.frame(marriage_period_start=c(1990,2000,2004),marriage_period_end=c(1995,2004,2008),
  ever_cohabited_before_first_marriage=c(.630,.745,.801),
  cohabited_with_spouse_before_first_marriage=c(.613,.725,.787),
  mean_premarital_cohabitation_years=c(2.7,3.7,3.8))
save_table(p,"great_britain_premarital_cohabitation.csv")

paths <- c(
  "/peoplepopulationandcommunity/birthsdeathsandmarriages/marriagecohabitationandcivilpartnerships/datasets/marriagesinenglandandwales2013/2019/marriages2019workbook.xlsx",
  "/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/backseriespopulationestimatesbymaritalstatusselectedyearsfrommid1971tomid2001singleyearofageandsexforenglandandwales/current/historicmaritalstatusestimatesmid1971tomid2001.xls",
  "/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesbymaritalstatusandlivingarrangements/current/previous/v7/2020ewpublicationtablefinal.xlsx",
  "/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/datasets/divorcesinenglandandwales/2019/divorcesworkbook2019fortracker.xlsx",
  "/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/bulletins/divorcesinenglandandwales/2014-02-06/60112210.xls",
  "/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/bulletins/divorcesinenglandandwales/2014-02-06/63abbfc5")
raw_files <- c("ons_marriages_2019.xlsx","ons_marital_population_1971_2001.xls",
  "ons_marital_population_2002_2019.xlsx","ons_divorces_2019.xlsx",
  "ons_divorce_probability_2010.xls","ons_divorce_death_cumulative_2010.csv")
urls <- paste0("https://www.ons.gov.uk/file?uri=",vapply(paths,URLencode,"",reserved=TRUE))
urls[6] <- paste0("https://www.ons.gov.uk/generator?format=csv&uri=",URLencode(paths[6],reserved=TRUE))
save_table(data.frame(file=raw_files,url=urls,retrieved="2026-11-02",
  licence="Open Government Licence v3.0",
  md5=unname(tools::md5sum(file.path(root,raw_files)))),"sources.csv")
message("Rebuilt England/Wales observations. No interpolation or scenario hazards selected.")
source(file.path(root,"prepare_marital_divorce_anchors.R"),local=TRUE)
