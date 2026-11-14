# rate messages keep their layout

    Code
      mortality(rates)
    Condition
      Error in `mortality()`:
      ! Age bands in `rates` must not overlap.
      i In calendar period [2000, Inf).
      x Bands [0, 20) and [10, Inf) overlap for sex "m".

---

    Code
      run_simulation(model, 2)
    Condition
      Error in `run_simulation()`:
      ! Mortality rates must cover the whole simulation.
      x Nothing applies at time 1.
      i The simulation runs from 0 to 2. Extend the time_start and time_end intervals to cover 1, or shorten `years`.

# union formation messages keep their layout

    Code
      union_formation("rate", rates)
    Condition
      Error in `union_formation()`:
      ! Column previous_union of `rates` must use pooled or detailed categories.
      x Categories: "none" and "separated".
      i Use "none" and "dissolved", or "none", "separated" and "widowed", in every group and period, without mixing the two.

