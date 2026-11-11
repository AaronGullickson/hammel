# Build a minimal population around union records for dissolution tests. Each
# union gets its own wife and husband unless `wife_id`/`husband_id` are given.
# Partners alive at -30 die at the union end for death endings. Optional
# `group_w`/`group_h` columns become the partners' declared `group`.
union_population <- function(unions, end_time) {
  n <- nrow(unions)
  if (!"wife_id" %in% names(unions)) {
    unions$wife_id <- 2L * seq_len(n) - 1L
  }
  if (!"husband_id" %in% names(unions)) {
    unions$husband_id <- 2L * seq_len(n)
  }
  ids <- sort(unique(c(unions$wife_id, unions$husband_id)))
  persons <- data.frame(
    person_id = ids,
    sex = c("m", "f")[1L + (ids %in% unions$wife_id)],
    birth_time = rep(-30, length(ids)),
    death_time = rep(NA_real_, length(ids)),
    mother_id = rep(NA_integer_, length(ids)),
    father_id = rep(NA_integer_, length(ids))
  )
  for (side in c("wife", "husband")) {
    dies <- unions$end_reason %in% paste0("death_", side)
    who <- match(unions[[paste0(side, "_id")]][dies], ids)
    persons$death_time[who] <- unions$end_time[dies]
  }
  groups <- NULL
  if ("group_w" %in% names(unions)) {
    persons$group <- rep(NA_character_, length(ids))
    persons$group[match(unions$wife_id, ids)] <- unions$group_w
    persons$group[match(unions$husband_id, ids)] <- unions$group_h
    unions$group_w <- NULL
    unions$group_h <- NULL
    groups <- "group"
  }
  hammel_population(persons, unions, end_time = end_time, groups = groups)
}
