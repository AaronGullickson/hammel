library(here)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(hammel)

pop <- read_fwf(
  here("data-raw", "pop_usa1880", "usa_00159.dat.gz"),
  col_positions = fwf_cols(
    SEX = c(85, 85),
    AGE = c(86, 88),
    RACE = c(90, 90),
    BPL = c(94, 96),
    MBPL = c(102, 104),
    FBPL = c(110, 112),
    MOMLOC = c(72, 73),
    POPLOC = c(74, 75),
    STEPPOP = c(76, 76),
    SPLOC = c(77, 78),
    RELATE = c(79, 80),
    SERIAL = c(11, 18),
    PERNUM = c(55, 58)
  ), 
  col_types = cols(.default = "i")
)

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

pop <- pop |>
  mutate(
    person_id = 1:nrow(pop), 
    sex = if_else(SEX == 1, "m", "f"),
    birth_time = 1880.5 - AGE + runif(n()),
    death_time = NA_real_,
    race = case_when(
      RACE == 1 ~ "White",
      RACE == 2 | RACE == 3 ~ "Black"
    ),
    birthplace = code_bpl(BPL),
    birthplace_mom = code_bpl(MBPL),
    birthplace_dad = code_bpl(FBPL),
    generation = case_when(
      birthplace != "native" ~ "1",
      birthplace_mom == "native" & birthplace_dad == "native" ~ "3+",
      birthplace_mom == "native" | birthplace_dad == "native" ~ "2.5",
      TRUE ~ "2"
    )
  )

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

# should be -1
cor(pop$AGE, pop$birth_time)

pop <- pop |>
  filter(!is.na(race), !is.na(birthplace), !is.na(generation))

mom_lookup <- pop |>
  select(SERIAL, PERNUM, person_id) |>
  rename(mother_id = person_id, MOMLOC = PERNUM)

dad_lookup <- pop |>
  select(SERIAL, PERNUM, person_id) |>
  rename(father_id = person_id, POPLOC = PERNUM)

pop <- pop |>
  left_join(mom_lookup, by = c("SERIAL", "MOMLOC")) |>
  left_join(dad_lookup, by = c("SERIAL", "POPLOC"))

pop <- pop |>
  select(person_id, sex, birth_time, death_time, mother_id, father_id,
         race, birthplace, generation)

pop_usa1880 <- list(persons = pop, unions = NULL, end_time = 1880.5)
# assign groups
pop_usa1880$groups <- list(
  race = hammel::group(
    categories = c("White", "Black"),
    inheritance = group_inheritance("hierarchy", order = c("Black", "White"))
  ),
  generation = hammel::group(
    categories = c("1", "1.5", "2", "2.5", "3+"),
    inheritance = group_inheritance("mother")
  )
)

# TODO: create arrival_time (and assign as 1.5)
# TODO: sample households
# TODO create custom generation inheritance
# TODO distribute by religion

pop_usa1880 <- do.call(hammel::hammel_population, pop_usa1880)
summary(pop_usa1880)
