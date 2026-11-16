# checks format their messages consistently

    Code
      set_rates(tibble::tibble(age_start = 0, rate = 1))
    Condition
      Error in `set_rates()`:
      ! `rates` is missing required columns.
      x Missing: age_end.

---

    Code
      set_rates(tibble::tibble(age_start = 0, age_end = 1, rate = c(1, NA)))
    Condition
      Error in `set_rates()`:
      ! Column rate of `rates` must contain finite numbers.
      x Missing at row 2.

---

    Code
      check_time("2000", "end_time")
    Condition
      Error:
      ! `end_time` must be a single finite number of years.
      x It is `"2000"`.

