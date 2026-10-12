library(hammel)

rates <- data.frame(age_start = c(0, 1), age_end = c(1, Inf),
                    rate = c(0.12, 0.24))
wait <- hammel:::hazard_wait_cpp(rates, 10/12, 0.05)
stopifnot(isTRUE(all.equal(as.numeric(wait), 3.5/12)))
constant <- data.frame(age_start = 0, age_end = Inf, rate = 0.2)
stopifnot(isTRUE(all.equal(as.numeric(hammel:::hazard_wait_cpp(
  constant, c(0, 10), c(0.2, 0.4))), c(1, 2))))
zero <- rates; zero$rate <- 0
stopifnot(is.infinite(hammel:::hazard_wait_cpp(zero, 0, 0.5)))
zero$rate[2] <- 1
stopifnot(hammel:::hazard_wait_cpp(zero, 0, 0.5) == 1.5)

# Survival agrees with the integrated hazard across the age boundary.
set.seed(42)
n <- 100000
waits <- hammel:::hazard_wait_cpp(rates, rep(0, n), rexp(n))
expected <- exp(-0.12 - 0.24)
stopifnot(abs(mean(waits > 2) - expected) <
            6 * sqrt(expected * (1 - expected) / n))
