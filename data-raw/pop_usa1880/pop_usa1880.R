# This script will generate a hammel_population sample of the 1880 Census sample
# from IPUMS. The script includes only the white and black population and
# collapses immigrant origins into NW Europe (including Germany), Irish,
# South Europe, East Europe, and Other. Arrival times for immigrants are
# sampled from a uniform distribution from the lastest non-native birth
# time for the family unit. Union start times are sampled from a uniform
# distribution from when the younger member was 15 or first birth time
# for the union.
#
# These data are generated for demonstration purposes only and should not be
# directly used for research purposes.

# TODO: create custom generation inheritance

# Loab libraries --------------------------------------------------------------

library(here)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(hammel)

set.seed(42)

# Create custom functions -----------------------------------------------------

# function for coding birthplace from IPUMS categories
code_bpl <- function(bpl) {

  case_when(
    bpl >= 900 ~ NA_character_,
    bpl < 100 ~ "native",
    bpl == 414 ~ "irish",
    (bpl >= 400 & bpl <= 429) | bpl == 453 ~ "nw_europe",
    bpl >= 430 & bpl <= 440 ~ "south_europe",
    bpl >= 450 & bpl <= 465 ~ "east_europe",
    TRUE ~ "other"
  )
}

# custom function for inheritance. We wrap it in local so it and its helper
# functions will travel with the saved object.
gen_inheritance <- local({
  # Report generation as numeric, rounded down.
  flatten_gen <- function(gen) {
    floor(as.numeric(sub("\\+$", "", gen)))
  }

  # Increment the generation and cap it at "3+".
  increment_gen <- function(gen) {
    gen <- min(3, gen + 1)
    if (gen == 3) "3+" else as.character(gen)
  }

  function(mother_group, father_group, ...) {
    gen_mother <- flatten_gen(mother_group)

    if (is.na(father_group) || mother_group == father_group) {
      return(increment_gen(gen_mother))
    }

    gen_father <- flatten_gen(father_group)

    # Exactly one first-generation parent produces generation 2.5.
    if (xor(gen_mother == 1, gen_father == 1)) {
      return("2.5")
    }

    increment_gen(min(gen_mother, gen_father))
  }
})

# Read in data and sample households ------------------------------------------

pop <- read_fwf(
  here("data-raw", "pop_usa1880", "usa_00161.dat.gz"),
  col_positions = fwf_cols(
    SEX = c(85, 85),
    AGE = c(86, 88),
    RACE = c(90, 90),
    BPL = c(94, 96),
    MBPL = c(102, 104),
    FBPL = c(110, 112),
    MOMLOC = c(71, 72),
    STEPMOM = c(73, 73),
    POPLOC = c(74, 75),
    STEPPOP = c(76, 76),
    SPLOC = c(77, 78),
    RELATE = c(79, 80),
    SERIAL = c(11, 18),
    PERNUM = c(55, 58),
    FAMUNIT = c(69, 70)
  ), 
  col_types = cols(.default = "i")
) |>
  # sample 10K households
  filter(SERIAL %in% sample(unique(SERIAL), size = 10000))

# Code basic individual demographic information -------------------------------

pop <- pop |>
  mutate(
    person_id = 1:nrow(pop), 
    sex = if_else(SEX == 1, "m", "f"),
    birth_time = 1880.5 - AGE - runif(n()),
    death_time = NA_real_,
    race = case_when(
      RACE == 1 ~ "White",
      RACE == 2 | RACE == 3 ~ "Black"
    ),
    birthplace = code_bpl(BPL),
    birthplace_mom = code_bpl(MBPL),
    birthplace_dad = code_bpl(FBPL),
    mom_native = birthplace_mom == "native",
    dad_native = birthplace_dad == "native",
    generation = case_when(
      birthplace != "native" ~ "1",
      coalesce(mom_native, dad_native) &
        coalesce(dad_native, mom_native) ~ "3+",
      coalesce(mom_native, dad_native) |
        coalesce(dad_native, mom_native) ~ "2.5",
      !is.na(mom_native) | !is.na(dad_native) ~ "2",
      TRUE ~ NA_character_
    )
  )

# checks
table(pop$SEX, pop$sex, exclude = NULL)
table(pop$RACE, pop$race, exclude = NULL)
table(pop$BPL, pop$birthplace, exclude = NULL)
table(pop$birthplace, pop$generation)
tf <- pop$birthplace == "native"
table(
  pop$birthplace_mom[tf],
  pop$birthplace_dad[tf],
  pop$generation[tf]
)

# should be close to -1
cor(pop$AGE, pop$birth_time)

# remove missing values
pop <- pop |>
  filter(!is.na(race), !is.na(birthplace), !is.na(generation))

# Add in mother, father, and spouse links -------------------------------------

# remove step parent links
pop <- pop |>
  mutate(
    POPLOC = if_else(STEPPOP > 0, NA, POPLOC),
    MOMLOC = if_else(STEPMOM > 0, NA, MOMLOC)
  )

mom_lookup <- pop |>
  select(SERIAL, PERNUM, person_id) |>
  rename(mother_id = person_id, MOMLOC = PERNUM)

dad_lookup <- pop |>
  select(SERIAL, PERNUM, person_id) |>
  rename(father_id = person_id, POPLOC = PERNUM)

partner_lookup <- pop |>
  select(SERIAL, PERNUM, person_id) |>
  rename(spouse_id = person_id, SPLOC = PERNUM)

pop <- pop |>
  left_join(mom_lookup, by = c("SERIAL", "MOMLOC")) |>
  left_join(dad_lookup, by = c("SERIAL", "POPLOC")) |>
  left_join(partner_lookup, by = c("SERIAL", "SPLOC"))

# Determine arrival time for non-natives --------------------------------------

# determine arrival time by runif between last non-native
# birth time and first-native birth/`now
immigrant_ids <- pop$person_id[which(pop$birthplace != "native")]

family_arrivals <- pop |>
  mutate(
    us_born_child = birthplace == "native" &
      (mother_id %in% immigrant_ids | father_id %in% immigrant_ids)
  ) |>
  group_by(SERIAL, FAMUNIT) |>
  summarize(
    last_foreign_birth = max(
      c(-Inf, birth_time[which(birthplace != "native")])
    ),
    first_native_birth = min(c(Inf, birth_time[which(us_born_child)])),
    .groups = "drop"
  ) |> 
  mutate(
    earliest_arrival = if_else(
      is.finite(last_foreign_birth), last_foreign_birth, NA_real_
    ),
    latest_arrival = if_else(
      is.finite(first_native_birth), first_native_birth, 1880.5
    ),
    # possible that some foreign born births come after native
    # so assign at exactly last foreign birth for these cases
    latest_arrival = pmax(earliest_arrival, latest_arrival),
    arrival_time = earliest_arrival +
      runif(n()) * (latest_arrival - earliest_arrival)
  ) |>
  select(SERIAL, FAMUNIT, arrival_time)

pop <- pop |>
  left_join(family_arrivals, by =c("SERIAL", "FAMUNIT")) |>
  mutate(
    # don't use arrival time for native born members of immigrant families
    arrival_time = if_else(birthplace == "native", NA, arrival_time),
    # check 1.5 generation status
    generation = if_else(
      generation == "1" & (arrival_time - birth_time) < 18, "1.5", generation
    )
  )

pop |>
  group_by(generation) |>
  summarize(arrival_time = mean(arrival_time))

# Create the unions dataset ---------------------------------------------------

# create  unions
unions <- pop |>
  filter(sex == "f") |>
  select(
    person_id,
    birth_time_w = birth_time,
    arrival_time_w = arrival_time,
    birthplace_w = birthplace,
  ) |>
  inner_join(
    pop |>
      filter(sex == "m") |>
      select(
        spouse_id,
        birth_time_h = birth_time,
        arrival_time_h = arrival_time,
        birthplace_h = birthplace,
        husband_id = person_id
      ),
    by = join_by(person_id == spouse_id)
  ) |>
  mutate(
    union_id = 1:n(),
    younger_birth = pmax(birth_time_w, birth_time_h),
    mixed_nativity = (birthplace_w == "native") != (birthplace_h == "native"),
    arrival_time = pmax(arrival_time_w, arrival_time_h, na.rm = TRUE)
  ) |>
  select(union_id, wife_id = person_id, husband_id, husband_id,
         younger_birth, mixed_nativity, arrival_time, starts_with("birthplace"))

# use birthdate of youngest spouse, first child birthdate from union, and arrival time
# to get range of potential starting times and then draw a starting time for each couple
unions <- unions |>
  left_join(
    pop |>
      select(mother_id, father_id, birth_time, birthplace) |>
      rename(birth_time_c = birth_time, birthplace_c = birthplace),
    by = join_by(wife_id == mother_id, husband_id == father_id)
  ) |>
  group_by(union_id, younger_birth, mixed_nativity, arrival_time) |>
  summarize(
    first_child_time = min(birth_time_c),
    first_child_born_abroad = any(
      birth_time_c == first_child_time & birthplace_c != "native",
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    early_date = pmin(younger_birth + 15, 1880),
    # check for children born before early date and modify
    early_date = if_else(
      !is.na(first_child_time) & first_child_time <= early_date,
      first_child_time - 1,
      early_date
    ),
    early_date = pmax(early_date, younger_birth),
    # check for mixed nativity - allow pre-migration marriage only if children
    # born abroad
    early_date = if_else(
      mixed_nativity &
        !(first_child_born_abroad & first_child_time < arrival_time),
      pmax(early_date, arrival_time),
      early_date
    ),
    late_date = coalesce(first_child_time, 1880.5),
    early_date = if_else(early_date <= late_date, early_date, NA_real_),
    start_time = early_date + runif(n()) * (late_date - early_date)
  ) |>
  select(union_id, start_time, arrival_time) |>
  right_join(unions) |>
  mutate(
    end_time = NA_real_, 
    end_reason = NA_character_,
    pre_migration = !is.na(arrival_time) & arrival_time > start_time
  ) |>
  select(union_id, wife_id, husband_id, start_time, end_time, end_reason, 
         pre_migration)

# Convert into hammel_population object ---------------------------------------

pop <- pop |>
  select(person_id, sex, birth_time, death_time, mother_id, father_id,
         birthplace, arrival_time, generation, race)

pop_usa1880 <- list(persons = pop, unions = unions, end_time = 1880.5)
# assign groups
pop_usa1880$groups <- list(
  race = hammel::group(
    categories = c("White", "Black"),
    inheritance = group_inheritance("hierarchy", order = c("Black", "White"))
  ),
  generation = hammel::group(
    categories = c("1", "1.5", "2", "2.5", "3+"),
    inheritance = group_inheritance("custom", fun = gen_inheritance)
  )
)

pop_usa1880 <- do.call(hammel::hammel_population, pop_usa1880)
summary(pop_usa1880)

save(pop_usa1880, file = here("data", "pop_usa1880.rda"))
