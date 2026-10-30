make_group_history <- function() {
  # The first child deliberately differs from its mother: fertility denominators
  # describe mothers, while population counts describe the individuals themselves.
  persons <- data.frame(person_id = c(10, 20, 30, 40, 50, 60),
    sex = c("f", "m", "f", "m", "m", "f"),
    birth_time = c(-20, -30, 0.5, 0.75, 0, -40),
    death_time = c(NA, 0.75, NA, NA, 0.5, 1),
    mother_id = c(NA, NA, 10, 10, NA, NA),
    father_id = c(NA, NA, 20, 20, NA, NA),
    group = c("A", "B", "B", "A", "B", "B"))
  unions <- data.frame(union_id = 1, wife_id = 10, husband_id = 20,
    start_time = 0.25, end_time = 0.75, end_reason = "death_husband")
  hammel_population(persons, unions, end_time = 1)
}

test_that("group Lexis counts preserve exact maternal and spouse classification", {
  p <- make_group_history()
  counts <- collect_demographic_counts(p, c(0, 1), c(0, Inf),
                                   by_union_type = TRUE, by_group = TRUE)
  expect_s3_class(counts, "tbl_df")
  expect_identical(unique(counts$group), c("A", "B"))
  a <- subset(counts, group == "A" & sex == "f")
  expect_equal(a$person_years, c(0.5, 0.5))
  expect_equal(a$births, c(1, 1))
  expect_equal(a$rate_fertility, c(2, 2))
  expect_false(any(c("unions", "rate_union") %in% names(a)))
  b <- subset(counts, group == "B" & sex == "m")
  expect_equal(b$person_years, c(0.75, 0.5))
  expect_equal(b$deaths, c(1, 1))
  expect_true(all(is.na(b$births)))
  expect_equal(subset(counts, group == "B" & sex == "f")$births, c(0, 0))
  unknown <- attr(counts, "unassigned_births")
  expect_s3_class(unknown, "tbl_df")
  expect_named(unknown, c("time_start", "time_end", "births"))
  expect_equal(unknown$births, 1)
})

test_that("group exposure and events conserve pooled cells across boundaries", {
  p <- make_group_history()
  for (by_union_type in c(FALSE, TRUE)) {
    counts <- collect_demographic_counts(p, c(0, 0.25, 0.5, 0.75, 1),
      c(0, 0.5, 20, 20.5, 30, 40, Inf), by_union_type, by_group = TRUE)
    pooled <- collect_demographic_counts(p, c(0, 0.25, 0.5, 0.75, 1),
      c(0, 0.5, 20, 20.5, 30, 40, Inf), by_union_type)
    a <- subset(counts, group == "A")
    b <- subset(counts, group == "B")
    for (column in intersect(c("person_years", "births", "deaths", "unions"), names(pooled)))
      expect_equal(a[[column]] + b[[column]], pooled[[column]])
    expect_identical(attr(counts, "unassigned_births"), attr(pooled, "unassigned_births"))
    expect_equal(a$rate_mortality[a$person_years > 0],
                 a$deaths[a$person_years > 0] / a$person_years[a$person_years > 0])
    expect_true(all(is.na(counts$rate_mortality[counts$person_years == 0])))
    expect_identical(collect_demographic_counts(p, c(0, 0.25, 0.5, 0.75, 1),
      c(0, 0.5, 20, 20.5, 30, 40, Inf), by_union_type, by_group = FALSE), pooled)
    if (by_union_type) {
      overall <- collect_demographic_counts(p, c(0, 0.25, 0.5, 0.75, 1),
        c(0, 0.5, 20, 20.5, 30, 40, Inf), by_group = TRUE)
      partnered <- subset(counts, union_type == "partnered")
      unpartnered <- subset(counts, union_type == "unpartnered")
      for (column in c("person_years", "births", "deaths"))
        expect_equal(partnered[[column]] + unpartnered[[column]], overall[[column]])
      expect_equal(subset(counts, group == "A" & union_type == "partnered" & births > 0)$time_start, 0.5)
      expect_equal(subset(counts, group == "A" & union_type == "unpartnered" & births > 0)$time_start, 0.75)
    }
  }
})

test_that("population sizes preserve time order and count each person's group", {
  p <- make_group_history()
  times <- c(1, 0, 0.5, 0.75, 0.5)
  sizes <- calculate_population_size(p, times, by_group = TRUE)
  expect_s3_class(sizes, "tbl_df")
  expect_identical(sizes$time, rep(times, each = 2))
  expect_identical(sizes$group, rep(c("A", "B"), length(times)))
  expect_equal(subset(sizes, group == "A")$population, c(2, 1, 1, 2, 1))
  expect_equal(subset(sizes, group == "B")$population, c(1, 3, 3, 2, 3))
  pooled <- calculate_population_size(p, times)
  expect_equal(rowSums(matrix(sizes$population, ncol = 2, byrow = TRUE)), pooled$population)
  expect_identical(calculate_population_size(p$persons, times, TRUE), sizes)
  expect_identical(calculate_population_size(p, times, FALSE), pooled)
})

test_that("numeric group types persist and bare histories keep maternal links", {
  p <- make_group_history()
  p$persons$group <- ifelse(p$persons$group == "A", 20L, 10L)
  counts <- collect_demographic_counts(p$persons, c(0, 1), c(0, Inf), by_group = TRUE)
  expect_type(counts$group, "integer")
  expect_identical(unique(counts$group), c(20L, 10L))
  expect_equal(subset(counts, group == 20L & sex == "f")$births, 2)
  sizes <- calculate_population_size(p, c(0, 1), TRUE)
  expect_type(sizes$group, "integer")
  expect_identical(sizes$group, c(20L, 10L, 20L, 10L))
})

test_that("TFR separates groups jointly with periods and union statuses", {
  grid <- expand.grid(age_start = c(15, 20), group = c("B", "A"),
    time_start = c(0, 1), union_type = c("overall", "partnered", "unpartnered"),
    stringsAsFactors = FALSE)
  grid$age_end <- grid$age_start + 5
  grid$time_end <- grid$time_start + 1
  grid$rate <- ifelse(grid$group == "A", 0.1, 0.2) *
    (grid$time_start + 1) * match(grid$union_type, c("overall", "partnered", "unpartnered"))
  result <- calculate_tfr(grid)
  expect_equal(nrow(result), 12)
  expect_equal(result$tfr, ifelse(result$group == "A", 1, 2) *
    (result$time_start + 1) * match(result$union_type, c("overall", "partnered", "unpartnered")))
  process_rates <- subset(grid, time_start == 0 & union_type != "overall")
  process_rates <- process_rates[c("age_start", "age_end", "group", "union_type", "rate")]
  names(process_rates)[names(process_rates) == "union_type"] <- "union_type"
  expect_equal(nrow(calculate_tfr(fertility(process_rates))), 4)
  grid$group <- ifelse(grid$group == "A", 20L, 10L)
  expect_type(calculate_tfr(grid)$group, "integer")
  # One group's invalid bands must fail even if another group fills the gap.
  bad <- grid
  bad$age_end[1] <- 19
  expect_error(calculate_tfr(bad), "contiguous")
  # Missing exposure affects only its own period/group/status estimate.
  incomplete <- grid
  incomplete$rate[1] <- NA_real_
  observed <- calculate_tfr(incomplete)
  affected <- observed$group == grid$group[1] &
    observed$time_start == grid$time_start[1] &
    observed$union_type == grid$union_type[1]
  expect_identical(is.na(observed$tfr), affected)
  expect_equal(observed$tfr[!affected], calculate_tfr(grid)$tfr[!affected])
})

test_that("group life tables follow analytical constant-hazard survival", {
  grid <- expand.grid(age_start = c(0, 5), group = c("B", "A"),
    sex = c("m", "f"), time_start = c(0, 1),
    union_type = c("overall", "partnered", "unpartnered"), stringsAsFactors = FALSE)
  grid$age_end <- ifelse(grid$age_start == 0, 5, Inf)
  grid$mx <- ifelse(grid$group == "A", 0.1, 0.2) *
    (grid$time_start + 1) * ifelse(grid$sex == "m", 2, 1) *
    match(grid$union_type, c("overall", "partnered", "unpartnered"))
  result <- calculate_life_table(grid, radix = 1000)
  expect_equal(nrow(result), nrow(grid))
  expect_equal(result$ex, 1/result$mx)
  expect_equal(result$lx[result$age_start == 0], rep(1000, 24))
  expect_equal(result$lx[result$age_start == 5],
               1000 * exp(-5 * result$mx[result$age_start == 5]))
  process_rates <- subset(grid, time_start == 0 & union_type == "overall")
  process_rates <- process_rates[c("age_start", "age_end", "group", "sex", "mx")]
  names(process_rates)[names(process_rates) == "mx"] <- "rate"
  expect_equal(nrow(calculate_life_table(mortality(process_rates))), 8)
  grid$group <- ifelse(grid$group == "A", 20L, 10L)
  expect_type(calculate_life_table(grid)$group, "integer")
})

test_that("Lexis group identifiers flow through derived measures", {
  counts <- collect_demographic_counts(make_group_history(), c(0, 1), c(0, 20, 40, Inf),
                                   by_union_type = TRUE, by_group = TRUE)
  tfr <- calculate_tfr(counts)
  expect_equal(nrow(tfr), 4)
  expect_identical(unique(tfr$group), c("A", "B"))
  expect_identical(attr(tfr, "unassigned_births"), attr(counts, "unassigned_births"))
  life <- calculate_life_table(counts)
  expect_equal(nrow(life), nrow(counts))
  expect_identical(unique(life$group), c("A", "B"))
})

test_that("group observations validate flags and required group data", {
  p <- make_group_history()
  for (bad in list(NULL, NA, 1, "yes", c(TRUE, FALSE), matrix(TRUE))) {
    expect_error(collect_demographic_counts(p, c(0, 1), c(0, Inf), by_group = bad), "by_group")
    expect_error(calculate_population_size(p, 0, by_group = bad), "by_group")
  }
  p$persons$group <- NULL
  expect_error(collect_demographic_counts(p, c(0, 1), c(0, Inf), by_group = TRUE), "group.*column")
  expect_error(calculate_population_size(p, 0, by_group = TRUE), "group.*column")
  expect_error(calculate_population_size(make_group_history(), 2, TRUE), "endpoint")
  expect_error(collect_demographic_counts(make_group_history(), c(0, 2), c(0, Inf), by_group = TRUE), "endpoint")
  rates <- data.frame(group = NA_character_, age_start = 0, age_end = Inf, rate = 0.1)
  expect_error(calculate_life_table(rates), "group")
  rates$age_end <- 5
  expect_error(calculate_tfr(rates), "group")
})

test_that("group observations return typed empty outputs", {
  p <- make_group_history()
  empty_times <- calculate_population_size(p, numeric(), by_group = TRUE)
  expect_s3_class(empty_times, "tbl_df")
  expect_equal(nrow(empty_times), 0)
  p$persons <- p$persons[FALSE, ]
  p$ancestry <- p$ancestry[FALSE, ]
  p$unions <- p$unions[FALSE, ]
  sizes <- calculate_population_size(p, c(0, 1), TRUE)
  expect_equal(nrow(sizes), 0)
  expect_type(sizes$group, "character")
  for (by_union_type in c(FALSE, TRUE)) {
    counts <- collect_demographic_counts(p, c(0, 1), c(0, Inf), by_union_type, TRUE)
    expect_s3_class(counts, "tbl_df")
    expect_equal(nrow(counts), 0)
    expect_type(counts$group, "character")
    expect_equal(attr(counts, "unassigned_births")$births, 0)
    expect_s3_class(calculate_tfr(counts), "tbl_df")
    expect_type(calculate_tfr(counts)$group, "character")
    expect_s3_class(calculate_life_table(counts), "tbl_df")
    expect_type(calculate_life_table(counts)$group, "character")
  }
})

test_that("groups with no survivors retain zero cells and boundary deaths", {
  p <- make_group_history()
  p$persons$death_time[p$persons$group == "B"] <- 0.75
  sizes <- calculate_population_size(p, c(0.75, 1), TRUE)
  expect_equal(subset(sizes, group == "B")$population, c(0L, 0L))
  counts <- collect_demographic_counts(p, c(0.75, 1), c(0, Inf), by_group = TRUE)
  b <- subset(counts, group == "B")
  expect_equal(b$person_years, c(0, 0))
  expect_equal(sum(b$deaths), 4)
  expect_true(all(is.na(b$rate_mortality)))
})
