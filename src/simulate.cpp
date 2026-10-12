#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

using namespace Rcpp;

// One piecewise constant age interval: [start, end), in years, with an annual
// person-year rate. An open final interval has end = infinity.
struct Band {
  double start, end, rate;
};

std::vector<Band> read_bands(DataFrame rates, const std::string& sex = "") {
  // Copy the R rate table into a growable C++ collection. Mortality selects
  // one sex; fertility uses all rows. Sort the copy without changing R inputs.
  NumericVector starts = rates["age_start"], ends = rates["age_end"];
  NumericVector values = rates["rate"];
  CharacterVector sexes;
  if (!sex.empty()) sexes = rates["sex"];
  std::vector<Band> bands;
  for (int i = 0; i < starts.size(); ++i) {
    if (sex.empty() || as<std::string>(sexes[i]) == sex)
      bands.push_back({starts[i], ends[i], values[i]});
  }
  std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) {
    return a.start < b.start;
  });
  return bands;
}

// Invert cumulative hazard through the known age schedule. A zero-rate band
// consumes time but no threshold. Absence of further hazard gives infinity.
double waiting_time(const std::vector<Band>& bands, double age, double threshold) {
  // threshold has an Exp(1) distribution, equivalent to -log(U). Consume it
  // using rate * elapsed years until the event lies within an age band.
  const double origin = age;
  for (const Band& band : bands) {
    if (band.end <= age) continue;
    age = std::max(age, band.start); // Skip uncovered ages without adding hazard.
    const double duration = band.end - age;
    if (band.rate > 0) {
      if (threshold <= band.rate * duration)
        // Include time crossed in previous bands, then the remaining time here.
        return age - origin + threshold / band.rate;
      threshold -= band.rate * duration;
    }
    age = band.end;
  }
  return std::numeric_limits<double>::infinity();
}

// Internal deterministic entry point for analytical validation of the sampler.
// [[Rcpp::export]]
NumericVector hazard_wait_cpp(DataFrame rates, NumericVector ages,
                             NumericVector thresholds) {
  if (ages.size() != thresholds.size()) stop("Sampler inputs must have equal length.");
  std::vector<Band> bands = read_bands(rates);
  NumericVector waits(ages.size());
  for (int i = 0; i < ages.size(); ++i) {
    if (!R_finite(ages[i]) || ages[i] < 0 || !R_finite(thresholds[i]) ||
        thresholds[i] <= 0) stop("Sampler ages and thresholds are invalid.");
    waits[i] = waiting_time(bands, ages[i], thresholds[i]);
  }
  return waits;
}

struct Person {
  // Public history fields plus working state. NA died means currently alive;
  // last_birth enforces spacing. male_slot locates this person in living males.
  double id, born, died, mother, father, last_birth;
  bool female;
  int male_slot;
};
// person is a zero-based position in people, not the public person_id.
// Each person retains only one next event; death=false means a birth.
struct Event { int person; bool death; };

double rate_integral(const std::vector<Band>& bands, double lo, double hi) {
  // Integrate the supplied annual rate over [lo, hi), accounting for bands
  // that only partly overlap the interval. Uncovered fertility ages add zero.
  double total = 0;
  for (const Band& band : bands) {
    if (band.end <= lo) continue;
    if (band.start >= hi) break;
    total += band.rate * std::max(0.0, std::min(hi, band.end) - std::max(lo, band.start));
  }
  return total;
}

std::vector<Band> adjusted_fertility(const std::vector<Band>& raw, double dt,
                                    int spacing, double max_age) {
  // With no spacing there is no eligibility correction. Otherwise build a
  // step-sized age schedule once, covering all ages reachable during this run.
  if (!spacing) return raw;
  double count = std::ceil(max_age / dt) + 1;
  if (count > std::numeric_limits<int>::max()) stop("Age grid is too large.");
  std::vector<double> probabilities;
  std::vector<Band> adjusted;
  double ineligible = 0;
  for (int k = 0; k < static_cast<int>(count); ++k) {
    if (k % 1024 == 0) checkUserInterrupt();
    const double lo = k * dt, hi = (k + 1.0) * dt;
    // Average the input rate over this age step. Estimate the fraction excluded
    // by recent births using the preceding spacing steps, then increase the
    // hazard among eligible women. This SOCSIM-style correction is approximate.
    const double rate = rate_integral(raw, lo, hi) / dt;
    const double eligible = 1 - ineligible;
    if (eligible <= 0) stop("Birth interval and fertility rates imply nonpositive eligibility.");
    adjusted.push_back({lo, hi, rate / eligible});
    // Numerically stable 1 - exp(-rate * dt). Add the new probability to the
    // moving window and remove the one leaving it before the next age step.
    probabilities.push_back(-std::expm1(-rate * dt));
    ineligible += probabilities.back();
    if (k >= spacing) ineligible -= probabilities[k - spacing];
  }
  return adjusted;
}

// [[Rcpp::export]]
List simulate_cpp(DataFrame persons, DataFrame mortality_rates,
                  DataFrame fertility_rates, double start, int steps,
                  double dt, int spacing) {
  // R performs input validation and seed handling. All engine state below is
  // local to this call; public times and rates remain in years and person-years.
  NumericVector ids = persons["person_id"], births = persons["birth_time"];
  NumericVector deaths = persons["death_time"], moms = persons["mother_id"];
  NumericVector dads = persons["father_id"];
  CharacterVector sexes = persons["sex"];
  // people retains every history, including deaths. males contains positions
  // of currently living males for fast uniform father selection. index maps
  // public IDs to positions while reconstructing initial maternal histories.
  std::vector<Person> people;
  std::vector<int> males;
  std::map<double, int> index;
  double next_id = 0, max_age = 0;
  for (int i = 0; i < ids.size(); ++i) {
    bool female = as<std::string>(sexes[i]) == "f";
    bool alive = NumericVector::is_na(deaths[i]);
    int slot = (!female && alive) ? static_cast<int>(males.size()) : -1;
    if (slot >= 0) males.push_back(i);
    people.push_back({ids[i], births[i], deaths[i], moms[i], dads[i],
                      -std::numeric_limits<double>::infinity(), female, slot});
    index[ids[i]] = i;
    next_id = std::max(next_id, ids[i]);
    if (alive) max_age = std::max(max_age, start - births[i]);
  }
  // On restart, recover the latest recorded maternal birth, including children
  // who have died. A mother with no recorded children has no spacing restriction.
  for (const Person& child : people) {
    if (!NumericVector::is_na(child.mother)) {
      auto found = index.find(child.mother);
      if (found != index.end())
        people[found->second].last_birth = std::max(people[found->second].last_birth, child.born);
    }
  }
  const auto male_rates = read_bands(mortality_rates, "m");
  const auto female_rates = read_bands(mortality_rates, "f");
  const auto fertility_rates_adjusted = adjusted_fertility(
    read_bands(fertility_rates), dt, spacing, max_age + steps * dt);
  // An ordered map stores only occupied steps, each with a growable event list.
  // It lets us jump over empty steps while retaining future events in memory.
  std::map<int, std::vector<Event>> queue;

  // These lambdas are small functions sharing the enclosing run's state.
  auto candidate = [&](double waiting, int step) -> int {
    // Floor the continuous waiting time to an integer step offset. The value
    // steps is a sentinel for no event before the half-open run endpoint.
    if (!std::isfinite(waiting)) return steps;
    double offset = std::floor(waiting / dt);
    if (offset >= steps - step) return steps;
    return step + static_cast<int>(offset);
  };
  auto threshold = [&]() {
    // Use R's RNG, shared with the seed policy in simulate(). Exp(1) is the
    // cumulative-hazard threshold; avoid zero to keep waiting times positive.
    double value;
    do { value = R::rexp(1.0); } while (value <= 0);
    return value;
  };
  auto schedule = [&](int i, int step) {
    const Person& p = people[i];
    if (!NumericVector::is_na(p.died)) return;
    double now = start + step * dt;
    const auto& death_rates = p.female ? female_rates : male_rates;
    // Sample independent death and birth candidates, incorporating future
    // age-band changes. Dead persons never acquire another pending event.
    int death_step = candidate(waiting_time(death_rates, now - p.born, threshold()), step);
    int birth_step = steps;
    if (p.female) {
      // Derive spacing from elapsed years; round only to resolve numerical
      // noise in already step-aligned birth dates, never calendar month counts.
      double remaining = p.last_birth + spacing * dt - now;
      double delay_steps = remaining > 0 ? std::ceil(remaining / dt - 1e-9) : 0;
      if (delay_steps < steps - step) {
        int first = step + static_cast<int>(delay_steps);
        double age = start + first * dt - p.born;
        // Start fertility sampling at the first step allowed by birth spacing.
        birth_step = candidate(waiting_time(fertility_rates_adjusted, age, threshold()), first);
      }
    }
    int due = std::min(death_step, birth_step);
    // Compare integer-step dates. Equality selects death, even if the continuous
    // birth candidate preceded death within that same step. Discard the loser.
    if (due < steps) queue[due].push_back({i, death_step <= birth_step});
  };
  // Initialize one next event for each living person. Known age changes are
  // already incorporated, so birthdays alone do not trigger resampling.
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    schedule(i, 0);
  }
  double birth_count = 0, death_count = 0;
  unsigned long long processed = 0;
  while (!queue.empty()) {
    // begin() is the earliest occupied step because std::map orders its keys.
    int step = queue.begin()->first;
    auto& events = queue.begin()->second;
    while (!events.empty()) {
      // Random sequential execution across people is distinct from the
      // death-before-birth tie rule for one person's competing candidates.
      // Drain this step fully, including new events scheduled into it.
      if (++processed % 1024 == 0) checkUserInterrupt();
      size_t choice = std::min(events.size() - 1,
        static_cast<size_t>(R::runif(0, 1) * events.size()));
      Event event = events[choice];
      // Remove without shifting the whole vector: fill the selected slot with
      // the last entry, then pop it. Event-list ordering has no meaning.
      events[choice] = events.back();
      events.pop_back();
      const int i = event.person;
      const double now = start + step * dt;
      if (event.death) {
        // Retain the person's history, record death, and stop their scheduling.
        people[i].died = now;
        if (!people[i].female) {
          // Remove a dead male in constant time, updating the moved male's slot.
          int slot = people[i].male_slot;
          int moved = males.back();
          males[slot] = moved;
          people[moved].male_slot = slot;
          males.pop_back();
          people[i].male_slot = -1;
        }
        ++death_count;
      } else {
        // All living females face age-specific fertility risk, regardless of
        // unions. Sample from all currently living males without age or kinship
        // restrictions; preserve the birth with unknown father if none exist.
        double father = NA_REAL;
        if (!males.empty()) {
          size_t selected = std::min(males.size() - 1,
            static_cast<size_t>(R::runif(0, 1) * males.size()));
          father = people[males[selected]].id;
        }
        if (next_id >= 9007199254740991.0 ||
            people.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
          stop("Population exceeds supported ID or row limits.");
        // IDs increase beyond existing histories. Store newborn parentage and
        // equal-probability sex; no union is created by this engine.
        double mother = people[i].id;
        people[i].last_birth = now;
        bool female = R::runif(0, 1) < 0.5;
        int child = static_cast<int>(people.size());
        int slot = female ? -1 : static_cast<int>(males.size());
        people.push_back({++next_id, now, NA_REAL, mother, father,
                          -std::numeric_limits<double>::infinity(), female, slot});
        if (!female) males.push_back(child);
        ++birth_count;
        // Both newborn and mother need fresh candidates. The mother now has a
        // new last_birth; spacing delays her next birth, but not a possible death.
        schedule(child, step);
        schedule(i, step);
      }
    }
    // All events in this step, including same-step follow-ups, have executed.
    queue.erase(step);
  }
  // Materialize ordinary R columns from the C++ records. The R wrapper restores
  // custom columns, closes existing unions affected by deaths, and assembles the
  // final population and run metadata. Pending events are not saved as a checkpoint.
  NumericVector out_id(people.size()), out_birth(people.size()), out_death(people.size());
  NumericVector out_mother(people.size()), out_father(people.size());
  CharacterVector out_sex(people.size());
  for (size_t i = 0; i < people.size(); ++i) {
    const Person& p = people[i];
    out_id[i] = p.id; out_birth[i] = p.born; out_death[i] = p.died;
    out_mother[i] = p.mother; out_father[i] = p.father;
    out_sex[i] = p.female ? "f" : "m";
  }
  return List::create(_["persons"] = DataFrame::create(
    _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
    _["death_time"] = out_death, _["mother_id"] = out_mother, _["father_id"] = out_father),
    _["births"] = birth_count, _["deaths"] = death_count);
}
