# group and matching messages keep their layout

    Code
      set_groups(grouped_model(), race = group(categories = "A"))
    Condition
      Error in `set_groups()`:
      ! Group race has values missing from its categories.
      x Undeclared: "B".
      i Declared categories: "A".

---

    Code
      run_simulation(model, 0)
    Condition
      Error in `run_simulation()`:
      ! Matching `coefficients` must be named for each model column.
      i In the matching model.
      x Missing: "female_raceB".
      x Unknown: "female_raceC".
      i The model columns are "female_raceB".

---

    Code
      run_simulation(hammel_model(), 1)
    Condition
      Error in `run_simulation()`:
      ! The model can't be simulated yet.
      x Missing: A starting population, from `set_population()`.
      x Missing: A mortality process, from `set_mortality()`.
      x Missing: A fertility process, from `set_fertility()`.
      x Missing: A union formation process, from `set_union_formation()`.
      i See `summary(model)` for the model's current contents.

