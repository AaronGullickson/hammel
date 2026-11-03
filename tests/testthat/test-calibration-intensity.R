test_that("queued cell intensities equal explicit individual scans after every change", {
  born <- c(-19.875, -20.125, -24.875, -30.125)
  groups <- c(0L, 0L, 1L, 1L)
  base <- data.frame(age_start = c(15, 20, 20.05, 25, 30),
    age_end = c(20, 20.05, 24, 30, 40), rate = c(.1, 1, 0, .3, .05))
  search <- do.call(rbind, lapply(0:1, function(g)
    do.call(rbind, lapply(0:2, function(p)
      transform(base, .hammel_group = g,
        previous_union = c("none", "separated", "widowed")[p + 1],
        rate = rate * (1 + g) * (1 + p))))))
  fertility <- do.call(rbind, lapply(0:1, function(g)
    data.frame(age_start = c(15, 21, 25, 31), age_end = c(20.02, 25, 30.05, 45),
      rate = c(.2, .1, .05, .01) * (1 + g), .hammel_group = g)))
  ops <- data.frame(time = c(rep(0, 4), .125, .125, born[1] + 20.05, .2, .5, .7, .8,
                             .875, 1, 1.25, 5, 25, 25),
    person = c(0:3, 0, 0, 0, 1, 2, 1, 3, 0, 0, 0, 0, 0, 1),
    action = c(rep("activate", 4), "inspect", "ready", "inspect", "remove", "ready",
      "activate", "remove", "inspect", "activate", "inspect", "inspect", "remove", "remove"),
    previous = c(0, 1, 2, 0, rep(0, 5), 2, 0, 0, 1, rep(0, 4)),
    ready = c(-Inf, -Inf, .5, -Inf, -Inf, .875, -Inf, -Inf, 1.25, .7,
              -Inf, -Inf, 1, rep(-Inf, 4)))
  set.seed(729); rng <- .Random.seed
  actual <- inspect_intensity_tracker_cpp(born, groups, search, fertility,
    ops$time, as.integer(ops$person), ops$action, as.integer(ops$previous), ops$ready, 2L)
  expect_identical(.Random.seed, rng)
  active <- rep(FALSE, length(born)); previous <- rep(0, length(born)); ready <- rep(-Inf, length(born))
  hazard <- function(rates, i, now) {
    rows <- rates$.hammel_group == groups[i] &
      born[i] + rates$age_start <= now & now < born[i] + rates$age_end
    if ("previous_union" %in% names(rates))
      rows <- rows & rates$previous_union == c("none", "separated", "widowed")[previous[i] + 1]
    sum(rates$rate[rows])
  }
  for (k in seq_len(nrow(ops))) {
    i <- ops$person[k] + 1
    if (ops$action[k] == "activate") {
      active[i] <- TRUE; previous[i] <- ops$previous[k]; ready[i] <- ops$ready[k]
    } else if (ops$action[k] == "remove") active[i] <- FALSE
    else if (ops$action[k] == "ready" && active[i]) ready[i] <- ops$ready[k]
    expected <- matrix(0, actual$bins * 2L, 2)
    for (j in which(active)) {
      cell <- groups[j] * actual$bins + floor(born[j] / .25) - actual$first + 1
      expected[cell, 1] <- expected[cell, 1] + hazard(search, j, ops$time[k])
      if (ops$time[k] >= ready[j])
        expected[cell, 2] <- expected[cell, 2] + hazard(fertility, j, ops$time[k])
    }
    expect_equal(actual$intensities[[k]], expected, tolerance = 1e-12, info = paste("operation", k))
    expect_lte(actual$pending[k], sum(active))
  }
  expect_equal(actual$pending[length(actual$pending)], 0L)
  expect_true(all(tail(actual$intensities, 1)[[1]] == 0))
})

test_that("replacing pending boundaries cannot grow the queue", {
  rates <- data.frame(age_start = 15, age_end = 50, rate = .2)
  times <- c(0, seq(.001, .1, .001), 1, 2)
  actions <- c("activate", rep("ready", 100), "remove", "inspect")
  got <- inspect_intensity_tracker_cpp(-20, 0L, rates, rates, times,
    rep(0L, length(times)), actions, rep(0L, length(times)),
    c(-Inf, seq(.001, .1, .001) + .75, -Inf, -Inf), 1L)
  expect_true(all(got$pending <= 1))
  expect_equal(tail(got$pending, 2), c(0L, 0L))
  expect_true(all(tail(got$intensities, 1)[[1]] == 0))
})
