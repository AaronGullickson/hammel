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

