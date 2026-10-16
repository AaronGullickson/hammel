#include <Rcpp.h>
#include "kinship.h"
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

// Copy the R rate table into a growable C++ collection. Mortality selects
// one sex; fertility uses all rows. Sort the copy without changing R inputs.
std::vector<Band> read_bands(DataFrame rates, const std::string& sex = "",
                             const std::string& status = "") {
  NumericVector starts = rates["age_start"], ends = rates["age_end"];
  NumericVector values = rates["rate"];
  CharacterVector sexes;
  if (!sex.empty()) sexes = rates["sex"];
  CharacterVector statuses;
  if (!status.empty()) statuses = rates["union_status"];
  std::vector<Band> bands;
  for (int i = 0; i < starts.size(); ++i) {
    if ((sex.empty() || as<std::string>(sexes[i]) == sex) &&
        (status.empty() || as<std::string>(statuses[i]) == status))
      bands.push_back({starts[i], ends[i], values[i]});
  }
  std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) {
    return a.start < b.start;
  });
  return bands;
}

// Calculate waiting time until the next event. The threshold is drawn
// from Exp(1), equivalent in distribution to -log(U), with U uniform
// on (0, 1). Subtract rate * duration across successive age bands.
// In the band where the remaining threshold is reached, divide it
// by the rate to obtain the remaining waiting time.
double waiting_time(const std::vector<Band>& bands, double age, double threshold) {
  const double origin = age;
  for (const Band& band : bands) {
    // if we are older than the end of this band, skip it
    if (band.end <= age) continue;
    // if we are entering a band higher than current age, then age up
    age = std::max(age, band.start);
    // determine how long we will spend in this band
    const double duration = band.end - age;
    // A zero-rate band consumes time but no threshold.
    if (band.rate > 0) {
      if (threshold <= band.rate * duration)
        // The event will occur within this band. Return the total waiting time.
        // Include time crossed in previous bands, then the remaining time here.
        return age - origin + threshold / band.rate;
      //if the event doesn't happen consume hazard from this band
      threshold -= band.rate * duration;
    }
    age = band.end;
  }
  // if we got here, we are out of bands so the event never happens, return infinity
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
  // last_birth enforces spacing. male_slot locates this person in available males.
  double id, born, died, mother, father, last_birth;
  bool female;
  int male_slot;
};
// person is a zero-based position in people, not the public person_id.
// A version invalidates queued events after a change in union status.
enum class EventType { death, birth, search };
struct Event { int person; EventType type; unsigned long long version; };

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

// Calculate one constant eligible hazard at each original band's midpoint.
// Integrate the unadjusted rates across the preceding birth interval, including
// overlaps with earlier bands and gaps. Event sampling uses the ordinary Band
// representation and shared constant-hazard waiting_time() function.
std::vector<Band> adjusted_fertility(const std::vector<Band>& raw,
                                    double birth_interval) {
  for (const Band& band : raw) {
    if (!std::isfinite(band.end)) stop("Fertility age bands must have finite ends.");
  }
  if (birth_interval == 0) return raw;
  std::vector<Band> adjusted;
  for (const Band& band : raw) {
    checkUserInterrupt();
    const double midpoint = band.start + (band.end - band.start) / 2;
    const double eligible = 1 - rate_integral(raw,
      std::max(0.0, midpoint - birth_interval), midpoint);
    if (!(eligible > 0))
      stop("Birth interval and fertility rates imply nonpositive midpoint eligibility.");
    const double rate = band.rate / eligible;
    if (!std::isfinite(rate)) stop("Adjusted fertility rate cannot be represented.");
    adjusted.push_back({band.start, band.end, rate});
  }
  return adjusted;
}

// Internal deterministic entry point to inspect adjustment in analytical tests.
// [[Rcpp::export]]
DataFrame adjusted_fertility_cpp(DataFrame rates, double birth_interval) {
  if (!R_finite(birth_interval) || birth_interval < 0)
    stop("Birth interval must be finite and nonnegative.");
  const auto bands = adjusted_fertility(read_bands(rates), birth_interval);
  NumericVector starts(bands.size()), ends(bands.size()), values(bands.size());
  for (size_t i = 0; i < bands.size(); ++i) {
    starts[i] = bands[i].start;
    ends[i] = bands[i].end;
    values[i] = bands[i].rate;
  }
  return DataFrame::create(_["age_start"] = starts, _["age_end"] = ends,
                          _["rate"] = values);
}

// Partial Fisher-Yates using reusable positions, restored after each sample.
// This avoids allocating a hash node for each candidate at every search.
// Workspace positions are always identity between calls, including when the
// living-male pool grows or shrinks. They index males, not fixed person IDs.
std::vector<int> sample_male_positions(int count, int requested,
                                     std::vector<int>& workspace) {
  const int take = requested == 0 ? count : std::min(count, requested);
  std::vector<int> sample;
  sample.reserve(take);
  if (take == count) {
    for (int j = 0; j < count; ++j) sample.push_back(j);
    return sample;
  }
  while (workspace.size() < static_cast<size_t>(count))
    workspace.push_back(static_cast<int>(workspace.size()));
  std::vector<int> swaps;
  swaps.reserve(take);
  for (int j = 0; j < take; ++j) {
    if (j % 1024 == 0) checkUserInterrupt();
    const int selected = j + std::min(count - j - 1,
      static_cast<int>(R::runif(0, 1) * (count - j)));
    std::swap(workspace[j], workspace[selected]);
    sample.push_back(workspace[j]);
    swaps.push_back(selected);
  }
  // Undo in reverse so the next search starts with identity positions again.
  for (int j = take - 1; j >= 0; --j)
    std::swap(workspace[j], workspace[swaps[j]]);
  return sample;
}

// Evaluate candidate utilities without rounding age differences. Common linear
// and raw quadratic formulas use an exact native path; other formulas call the
// R model-matrix evaluator with candidate data only, never mutable engine state.
NumericVector matching_utilities(const NumericVector& gaps, const List& matching) {
  NumericVector quadratic = matching["quadratic"];
  NumericVector utility(gaps.size());
  if (quadratic.size() == 2) {
    const bool has_squared = as<bool>(matching["has_squared"]);
    for (R_xlen_t j = 0; j < gaps.size(); ++j) {
      const double squared = has_squared ? gaps[j] * gaps[j] : 0;
      if (!R_finite(squared)) stop("Matching formula must produce finite values for every candidate.");
      utility[j] = quadratic[0] * gaps[j] + quadratic[1] * squared;
    }
  } else {
    Function evaluate = matching["evaluate"];
    utility = as<NumericVector>(evaluate(gaps));
  }
  if (utility.size() != gaps.size()) stop("Matching formula returned the wrong number of utilities.");
  for (R_xlen_t j = 0; j < utility.size(); ++j) {
    if (!R_finite(gaps[j]) || !R_finite(utility[j]))
      stop("Matching ages and utilities must be finite.");
  }
  return utility;
}

NumericVector matching_probabilities(const NumericVector& utility) {
  if (utility.size() == 0) return NumericVector();
  const double maximum = *std::max_element(utility.begin(), utility.end());
  NumericVector probabilities(utility.size());
  double total = 0;
  for (R_xlen_t j = 0; j < utility.size(); ++j) {
    probabilities[j] = std::exp(utility[j] - maximum);
    total += probabilities[j];
  }
  for (double& probability : probabilities) probability /= total;
  return probabilities;
}

// One uniform draw selects a candidate from the softmax cumulative distribution.
int draw_partner_position(const NumericVector& probabilities) {
  if (probabilities.size() == 0) return -1;
  const double draw = R::runif(0, 1);
  double cumulative = 0;
  int last_positive = -1;
  for (R_xlen_t j = 0; j < probabilities.size(); ++j) {
    if (probabilities[j] > 0) last_positive = static_cast<int>(j);
    cumulative += probabilities[j];
    if (draw < cumulative) return static_cast<int>(j);
  }
  // Roundoff must not select a zero-weight tail candidate.
  return last_positive;
}

// Apply search eligibility after sampling; retain positions in the male pool.
// Further pair restrictions can be added here without changing sample size.
std::vector<int> screen_male_candidates(const std::vector<int>& males,
                                       const std::vector<int>& sample,
                                       const std::vector<Person>& people,
                                       double now, double min_male_age, int woman,
                                       hammel::Kinship& kinship) {
  std::vector<int> eligible;
  eligible.reserve(sample.size());
  const bool check_kinship = kinship.generations() > 0;
  for (int position : sample) {
    if (now - people[males[position]].born >= min_male_age &&
        (!check_kinship || !kinship.known_related(woman, males[position])))
      eligible.push_back(position);
  }
  return eligible;
}

// Select a partner from sampled positions in males, returning a people position.
// Age(male)-age(female) equals birth_time(female)-birth_time(male).
int select_partner(const std::vector<int>& males, const std::vector<int>& sample,
                   const std::vector<Person>& people, int woman, const List& matching) {
  if (sample.empty()) return -1;
  NumericVector gaps(sample.size());
  for (size_t j = 0; j < sample.size(); ++j)
    gaps[j] = people[woman].born - people[males[sample[j]]].born;
  const NumericVector utility = matching_utilities(gaps, matching);
  const NumericVector probabilities = matching_probabilities(utility);
  return males[sample[draw_partner_position(probabilities)]];
}

// Internal analytical/sampling diagnostics use the same scoring and draw helpers.
// [[Rcpp::export]]
List inspect_matching_cpp(NumericVector age_difference, List matching, int draws) {
  if (draws < 0) stop("Number of draws must be nonnegative.");
  const NumericVector utility = matching_utilities(age_difference, matching);
  const NumericVector probabilities = matching_probabilities(utility);
  IntegerVector choices(draws);
  for (int j = 0; j < draws; ++j) {
    if (j % 1024 == 0) checkUserInterrupt();
    const int selected = draw_partner_position(probabilities);
    choices[j] = selected < 0 ? NA_INTEGER : selected + 1;
  }
  return List::create(_["utilities"] = utility, _["probabilities"] = probabilities,
                      _["choices"] = choices);
}

// Internal sampling entry point for checking size, uniqueness, and inclusion.
// [[Rcpp::export]]
IntegerVector sample_male_positions_cpp(int count, int requested) {
  if (count < 0 || requested < 0) stop("Sample sizes must be nonnegative.");
  std::vector<int> workspace;
  return wrap(sample_male_positions(count, requested, workspace));
}

struct UnionRecord {
  double id, start, end;
  int wife, husband;
  bool wife_died = false;
};

// [[Rcpp::export]]
List simulate_cpp(DataFrame persons, DataFrame mortality_rates,
                  DataFrame fertility_rates, double start, double years,
                  double birth_interval, DataFrame unions,
                  int sample_size, List matching, double min_male_age,
                  int kinship_generations, bool rate_initiation,
                  DataFrame union_rates, List father_matching) {
  // R performs input validation and seed handling. All engine state below is
  // local to this call; public times and rates remain in years and person-years.
  NumericVector ids = persons["person_id"], births = persons["birth_time"];
  NumericVector deaths = persons["death_time"], moms = persons["mother_id"];
  NumericVector dads = persons["father_id"];
  CharacterVector sexes = persons["sex"];
  // people retains every history, including deaths. males contains positions
  // of living, unpartnered males for candidate sampling. index maps
  // public IDs to positions while reconstructing initial maternal histories.
  std::vector<Person> people;
  std::vector<int> males;
  std::vector<int> sampling_workspace;
  std::map<double, int> index;
  double next_id = 0;
  for (int i = 0; i < ids.size(); ++i) {
    bool female = as<std::string>(sexes[i]) == "f";
    bool alive = NumericVector::is_na(deaths[i]);
    int slot = (!female && alive) ? static_cast<int>(males.size()) : -1;
    if (slot >= 0) males.push_back(i);
    people.push_back({ids[i], births[i], deaths[i], moms[i], dads[i],
                      -std::numeric_limits<double>::infinity(), female, slot});
    index[ids[i]] = i;
    next_id = std::max(next_id, ids[i]);
  }
  std::vector<std::array<int, 2>> parent_links(people.size());
  for (int person = 0; person < static_cast<int>(people.size()); ++person) {
    parent_links[person] = {
      NumericVector::is_na(people[person].mother) ? -1 : index.at(people[person].mother),
      NumericVector::is_na(people[person].father) ? -1 : index.at(people[person].father)
    };
  }
  hammel::Kinship kinship(std::move(parent_links), kinship_generations);
  std::vector<UnionRecord> union_history;
  std::vector<std::vector<int>> person_unions(people.size());
  std::vector<int> current_partner(people.size(), -1);
  // Swap removal keeps pool updates constant-time and slots valid after moves.
  auto remove_male = [&](int man) {
    const int slot = people[man].male_slot;
    if (slot < 0) return; // Partnered men are already absent when they die.
    const int moved = males.back();
    males[slot] = moved;
    people[moved].male_slot = slot;
    males.pop_back();
    people[man].male_slot = -1;
  };
  auto add_male = [&](int man) {
    if (people[man].female || !NumericVector::is_na(people[man].died) ||
        current_partner[man] >= 0 || people[man].male_slot >= 0) return;
    people[man].male_slot = static_cast<int>(males.size());
    males.push_back(man);
  };
  double next_union_id = 0;
  NumericVector uids = unions["union_id"], wives = unions["wife_id"];
  NumericVector husbands = unions["husband_id"], starts = unions["start_time"];
  NumericVector ends = unions["end_time"];
  for (int u = 0; u < uids.size(); ++u) {
    const int woman = index.at(wives[u]), man = index.at(husbands[u]);
    union_history.push_back({uids[u], starts[u], ends[u], woman, man});
    next_union_id = std::max(next_union_id, uids[u]);
    if (!NumericVector::is_na(ends[u])) continue;
    if (!NumericVector::is_na(people[woman].died) ||
        !NumericVector::is_na(people[man].died) ||
        !people[woman].female || people[man].female)
      stop("Unions require living female/male partners in active input unions.");
    if (current_partner[woman] >= 0)
      stop("Union input cannot have multiple ongoing unions for a woman.");
    if (current_partner[man] >= 0)
      stop("Union input cannot have multiple ongoing unions for a man.");
    current_partner[woman] = man;
    current_partner[man] = woman;
    remove_male(man);
    person_unions[woman].push_back(u);
    person_unions[man].push_back(u);
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
  // All father and partner searches share sampling, eligibility and matching.
  // Only form_union changes partner state or removes the selected man.
  auto find_candidate = [&](int woman, double now, const List& choice_model) {
    if (males.empty()) return -1;
    const auto sample = sample_male_positions(static_cast<int>(males.size()), sample_size,
                                             sampling_workspace);
    const auto eligible = screen_male_candidates(males, sample, people, now, min_male_age,
                                                 woman, kinship);
    return select_partner(males, eligible, people, woman, choice_model);
  };
  auto form_union = [&](int woman, double now) {
    int partner = current_partner[woman];
    if (partner >= 0) return partner;
    partner = find_candidate(woman, now, matching);
    if (partner >= 0) {
      if (next_union_id >= 9007199254740991.0 ||
          union_history.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Union history exceeds supported ID or row limits.");
      const int u = static_cast<int>(union_history.size());
      union_history.push_back({++next_union_id, now, NA_REAL, woman, partner});
      person_unions[woman].push_back(u);
      person_unions[partner].push_back(u);
      current_partner[woman] = partner;
      current_partner[partner] = woman;
      remove_male(partner);
    }
    return partner;
  };
  const auto male_rates = read_bands(mortality_rates, "m");
  const auto female_rates = read_bands(mortality_rates, "f");
  const bool status_fertility = fertility_rates.containsElementNamed("union_status");
  const auto partnered_fertility = adjusted_fertility(
    read_bands(fertility_rates, "", status_fertility ? "partnered" : ""), birth_interval);
  const auto unpartnered_fertility = status_fertility ? adjusted_fertility(
    read_bands(fertility_rates, "", "unpartnered"), birth_interval) : partnered_fertility;
  const auto search_rates = rate_initiation ? read_bands(union_rates) : std::vector<Band>{};
  std::vector<unsigned long long> event_versions(people.size(), 0);
  // Continuous dates order events. Equivalent keys retain insertion order.
  std::multimap<double, Event> queue;
  const double end = start + years;
  auto candidate = [&](double waiting, double now) {
    double due = now + waiting;
    if (std::isfinite(waiting) && due <= now)
      stop("Event time cannot advance at this numerical scale.");
    return due;
  };
  auto threshold = [&]() {
    // Use R's RNG, shared with the seed policy in simulate(). Exp(1) is the
    // cumulative-hazard threshold; avoid zero to keep waiting times positive.
    double value;
    do { value = R::rexp(1.0); } while (value <= 0);
    return value;
  };
  auto schedule = [&](int i, double now) {
    const auto version = ++event_versions[i];
    const Person& p = people[i];
    if (!NumericVector::is_na(p.died)) return;
    const auto& death_rates = p.female ? female_rates : male_rates;
    double death_time = candidate(waiting_time(death_rates, now - p.born, threshold()), now);
    double birth_time = std::numeric_limits<double>::infinity();
    if (p.female) {
      const auto& fertility = current_partner[i] >= 0 ? partnered_fertility : unpartnered_fertility;
      double first = std::max(now, p.last_birth + birth_interval);
      if (first < end)
        birth_time = candidate(waiting_time(fertility,
          first - p.born, threshold()), first);
    }
    double search_time = std::numeric_limits<double>::infinity();
    if (rate_initiation && p.female && current_partner[i] < 0)
      search_time = candidate(waiting_time(search_rates, now - p.born, threshold()), now);
    const double due = std::min(death_time, std::min(birth_time, search_time));
    // Exact same-person ties select death, then birth, then search.
    const auto type = death_time <= due ? EventType::death :
      (birth_time <= search_time ? EventType::birth : EventType::search);
    if (due < end) queue.emplace(due, Event{i, type, version});
  };
  // Initialize one next event for each living person. Known age changes are
  // already incorporated, so birthdays alone do not trigger resampling.
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    schedule(i, start);
  }
  double birth_count = 0, death_count = 0;
  unsigned long long processed = 0;
  // Process events in chronological order until no scheduled events remain.
  while (!queue.empty()) {
    if (++processed % 1024 == 0) checkUserInterrupt();
    auto next = queue.begin();
    const double now = next->first;
    const Event event = next->second;
    queue.erase(next);
    const int i = event.person;
    if (event.version != event_versions[i]) continue;
    if (event.type == EventType::death) {
      // Retain the person's history, record death, and stop their scheduling.
      people[i].died = now;
      for (int u : person_unions[i]) {
        UnionRecord& relationship = union_history[u];
        if (!NumericVector::is_na(relationship.end)) continue;
        relationship.end = now;
        relationship.wife_died = i == relationship.wife;
        current_partner[relationship.wife] = -1;
        current_partner[relationship.husband] = -1;
        add_male(relationship.husband);
        // Her old birth candidate is invalid once her husband has died.
        if (rate_initiation && i == relationship.husband)
          schedule(relationship.wife, now);
      }
      if (!people[i].female) remove_male(i);
      ++death_count;
    } else if (event.type == EventType::search) {
      form_union(i, now);
      schedule(i, now);
    } else {
      double father = NA_REAL;
      int partner = current_partner[i];
      if (!rate_initiation && partner < 0) partner = form_union(i, now);
      // A nonmarital father stays available: parentage alone is not a union.
      const int father_position = partner >= 0 ? partner :
        (rate_initiation ? find_candidate(i, now, father_matching) : -1);
      if (father_position >= 0) father = people[father_position].id;
      // If no sampled men are eligible, preserve the birth with unknown father and no union.
      if (next_id >= 9007199254740991.0 ||
          people.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        stop("Population exceeds supported ID or row limits.");
      // IDs increase beyond existing histories. Store newborn parentage and
      // sex; union formation above happens before the newborn exists.
      double mother = people[i].id;
      people[i].last_birth = now;
      // Use 0.4886 as ffab - its the signature of Berkeley training!
      bool female = R::runif(0, 1) < 0.4886;
      int child = static_cast<int>(people.size());
      int slot = female ? -1 : static_cast<int>(males.size());
      people.push_back({++next_id, now, NA_REAL, mother, father,
                        -std::numeric_limits<double>::infinity(), female, slot});
      kinship.append(i, father_position);
      if (!female) males.push_back(child);
      person_unions.emplace_back();
      current_partner.push_back(-1);
      event_versions.push_back(0);
      ++birth_count;
      // Both newborn and mother need fresh candidates. The mother now has a
      // new last_birth; spacing delays her next birth, but not a possible death.
      schedule(child, now);
      schedule(i, now);
    }
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
  NumericVector out_uid(union_history.size()), out_wife(union_history.size());
  NumericVector out_husband(union_history.size()), out_start(union_history.size());
  NumericVector out_end(union_history.size());
  CharacterVector out_reason(union_history.size());
  CharacterVector input_reasons = unions["end_reason"];
  for (size_t u = 0; u < union_history.size(); ++u) {
    const UnionRecord& relationship = union_history[u];
    out_uid[u] = relationship.id;
    out_wife[u] = people[relationship.wife].id;
    out_husband[u] = people[relationship.husband].id;
    out_start[u] = relationship.start;
    out_end[u] = relationship.end;
    if (NumericVector::is_na(relationship.end)) out_reason[u] = NA_STRING;
    else if (u < static_cast<size_t>(unions.nrows())) {
      NumericVector input_ends = unions["end_time"];
      if (!NumericVector::is_na(input_ends[u])) out_reason[u] = input_reasons[u];
      else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
    } else out_reason[u] = relationship.wife_died ? "death_wife" : "death_husband";
  }
  return List::create(_["persons"] = DataFrame::create(
    _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
    _["death_time"] = out_death, _["mother_id"] = out_mother, _["father_id"] = out_father),
    _["unions"] = DataFrame::create(_["union_id"] = out_uid,
      _["wife_id"] = out_wife, _["husband_id"] = out_husband,
      _["start_time"] = out_start, _["end_time"] = out_end,
      _["end_reason"] = out_reason),
    _["births"] = birth_count, _["deaths"] = death_count);
}
