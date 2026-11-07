# Ballpark age-specific migration rates between US census regions, around
# 2000. Run from the package root. See usa_2000_regional_migration.README.md.
#
# Inputs, transcribed from U.S. Census Bureau, "Geographical Mobility: March
# 1999 to March 2000" (P20-538, Schachter 2001; the PDF is in this folder):
# - Figure 3: region-to-region migrants, thousands. The figure's values were
#   assigned to flows using Table D's regional in- and out-migrant totals,
#   which they reproduce exactly (checked below).
# - Table B: percent moving to a different state, by age (proxy for the age
#   pattern of moves between regions).
# Region populations are Census 2000 counts (thousands).

regions <- c("NE", "MW", "S", "W")
flows <- matrix(c(NA,  83, 405, 127,
                  72,  NA, 364, 204,
                  183, 416, NA, 432,
                  108, 223, 489, NA), 4, byrow = TRUE, dimnames = list(regions, regions))
stopifnot(rowSums(flows, na.rm = TRUE) == c(615, 640, 1031, 820),   # Table D outmigrants
          colSums(flows, na.rm = TRUE) == c(363, 722, 1258, 763))   # Table D inmigrants
population <- c(NE = 53594, MW = 64393, S = 100237, W = 63198)

# Table B, percent moving to a different state, by age; total 3.1. The first
# band (ages 1-4 in the source) is applied from age 0.
interstate <- data.frame(
  age_start = c(0, 5, 10, 20, 25, 30, 35, 45, 55, 65, 85),
  age_end   = c(5, 10, 20, 25, 30, 35, 45, 55, 65, 85, Inf),
  percent   = c(4.0, 3.3, 3.1, 6.2, 6.0, 4.0, 3.0, 1.9, 1.9, 1.1, 1.0))
multiplier <- interstate$percent / 3.1

pairs <- expand.grid(origin = regions, destination = regions, stringsAsFactors = FALSE)
pairs <- pairs[pairs$origin != pairs$destination, ]
pairs$all_ages <- flows[cbind(pairs$origin, pairs$destination)] / population[pairs$origin]
out <- merge(pairs, data.frame(interstate[c("age_start", "age_end")], multiplier = multiplier), by = NULL)
out$rate <- signif(out$all_ages * out$multiplier, 4)
out <- out[order(out$origin, out$destination, out$age_start), c("origin", "destination", "age_start", "age_end", "rate")]
write.csv(out, "validation/data/usa_2000_2001/usa_2000_regional_migration.csv", row.names = FALSE)
